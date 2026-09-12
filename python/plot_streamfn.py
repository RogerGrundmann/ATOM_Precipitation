# Plot the meridional mass streamfunction Psi(lat, z) for two iterations side by side
# from the write_meridional_streamfunction() CSV dumps, on a shared color scale, so the
# vertical (depth) collapse of the overturning is directly visible.
#
# ⚠ WHICH COLUMN. The CSV carries TWO streamfunctions and they disagree about closure by a
# factor of 50 to 500. `psi_kg_per_s` divides each level by n, the FLUID-cell count at that
# level, which varies with height wherever terrain varies with longitude, so Psi(ground) does
# not vanish even when every column's mass flux integral does. `psi_fixdiv_kg_per_s` integrates
# each column and then averages -- the definition -- and closes to 0.0000 on the closed-cell IC.
# The fixed column was added as an ADDITIONAL column so that no recorded Psi number changed
# meaning, and this script was never switched over: it plotted the broken one from the day the
# fix landed. At iteration 200 of output_pdc600 the plotted closure |Psi(gnd)|/cell read
# 0.095 / 0.135 / 0.049 / 0.074 / 0.008 / 0.462 at 75N/45N/15N/15S/45S/75S where the real
# figures are 0.0004 / 0.0008 / 0.0001 / 0.0009 / 0.0004 / 0.0048 -- i.e. the plot showed open
# cells, worst at the south pole, in a run whose cells are closed. Same shape as `u-v-Cell`
# versus `uv_plot` for the glyphs: the repaired quantity sat beside the broken one and the tool
# kept reading the broken one. DEFAULT IS NOW THE FIXED COLUMN, and the closure of BOTH is
# printed every run so this cannot recur silently.
#
# ⚠ WHICH SCALE. A SHARED LINEAR scale cannot show this model's weak cells at all. The six-cell
# structure spans a factor of 13 in amplitude -- at iteration 200 of output_pdc600 the Hadley
# cells are 97e9 kg/s and the polar ones 11.5e9 (75N) and 7.5e9 (75S) -- so with 21 levels over
# +-109 the contour interval is 10.9e9 and the polar cells span 1.06 and 0.69 INTERVALS: 75S gets
# not a single contour line drawn, and 85 deg spans 0.1-0.2. The cells are there and closed (sign
# reversals at 0, +-30, +-60 at every height; closure 0.0001-0.005 on psi_fixdiv); the plot simply
# could not resolve them. DEFAULT IS NOW SYMMETRIC-LOG levels, which render all six, and the
# Psi = 0 contour -- the cell BOUNDARY, i.e. what "closed" means visually -- is drawn heavier.
# --levels=lin restores the old linear spacing.
#
# usage: python3 plot_streamfn.py [outdir] [it1] [it2] [--col=fixed|old|both] [--levels=symlog|lin]
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, Normalize

COL_FIXED, COL_OLD = "psi_fixdiv_kg_per_s", "psi_kg_per_s"
argv  = [a for a in sys.argv[1:] if not a.startswith("--")]
opts  = [a for a in sys.argv[1:] if a.startswith("--")]
which = next((o.split("=", 1)[1] for o in opts if o.startswith("--col=")), "fixed")
scale = next((o.split("=", 1)[1] for o in opts if o.startswith("--levels=")), "symlog")

OUT   = argv[0] if argv else "output_0Ma"
iters = [int(a) for a in argv[1:3]] if len(argv) >= 3 else [20, 100]
BANDS = [75, 45, 15, -15, -45, -75]


def grid(df, col, lat, z):
    P = df.pivot(index="height_m", columns="lat_deg", values=col)
    return P.reindex(index=z, columns=lat).values / 1.0e9      # -> 1e9 kg/s


def load(it):
    df  = pd.read_csv(f"{OUT}/meridional_streamfunction_{it}.csv")
    lat = np.sort(df["lat_deg"].unique())          # ascending: 90S -> 90N
    z   = np.sort(df["height_m"].unique())         # ascending: surface -> lid
    if COL_FIXED not in df.columns:
        print(f"  ⚠ {COL_FIXED} absent (CSV predates the divisor fix) -- falling back to {COL_OLD}")
        return lat, z, {COL_OLD: grid(df, COL_OLD, lat, z)}
    return lat, z, {c: grid(df, c, lat, z) for c in (COL_FIXED, COL_OLD)}


def closure(lat, P):
    """|Psi(ground)| / max|Psi - Psi(ground)| per band -- 0 is a closed cell."""
    out = []
    for b in BANDS:
        j = int(np.argmin(np.abs(lat - b)))
        col = P[:, j]
        if not np.isfinite(col).any():
            out.append(float("nan")); continue
        g   = col[0]                                # z ascending, so row 0 is the ground
        dev = col - g
        cell = np.nanmax(np.abs(dev))
        out.append(abs(g) / cell if cell > 0 else float("nan"))
    return out


data = [load(it) for it in iters]
cols = [c for c in (COL_FIXED, COL_OLD) if c in data[0][2]]
plot_cols = cols if which == "both" else \
            [COL_OLD if which == "old" and COL_OLD in cols else cols[0]]

# The closure of BOTH columns, every run, so the two can never be confused again.
print(f"closure |Psi(ground)|/cell   (0 = closed)      bands: "
      + "  ".join(f"{b:>7}" for b in BANDS))
for it, (lat, z, P) in zip(iters, data):
    for c in cols:
        mark = " <- plotted" if c in plot_cols else ""
        print(f"  iter {it:<4} {c:<20} " + "  ".join(f"{v:7.4f}" for v in closure(lat, P[c])) + mark)

for col in plot_cols:
    vmax = max(np.nanmax(np.abs(P[col])) for _, _, P in data)
    if scale == "lin":
        levels, norm = np.linspace(-vmax, vmax, 21), Normalize(-vmax, vmax)
    else:
        # geometric spacing from lt to vmax: resolves a 13x amplitude range in one panel
        lt  = max(vmax / 300.0, 0.2)
        pos = np.geomspace(lt, vmax, 9)
        levels = np.concatenate([-pos[::-1], [0.0], pos])
        norm   = SymLogNorm(linthresh=lt, vmin=-vmax, vmax=vmax, base=10)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
    for ax, it, (lat, z, P) in zip(axes, iters, data):
        cf = ax.contourf(lat, z / 1000.0, P[col], levels=levels, norm=norm,
                         cmap="RdBu_r", extend="both")
        ax.contour(lat, z / 1000.0, P[col], levels=levels, colors="k", linewidths=0.3, alpha=0.5)
        # Psi = 0 is the CELL BOUNDARY -- the thing "closed cells" refers to.
        ax.contour(lat, z / 1000.0, P[col], levels=[0.0], colors="k", linewidths=1.3)
        ax.set_title(f"iter {it}")
        ax.set_xlabel("latitude [deg N]")
        ax.set_xlim(-90, 90)
        ax.set_xticks(np.arange(-90, 91, 30))
        ax.axvline(0, color="0.4", lw=0.6, ls="--")
    axes[0].set_ylabel("height [km]")
    cb = fig.colorbar(cf, ax=axes, shrink=0.9, pad=0.02)
    cb.set_label(r"$\Psi$  [$10^{9}$ kg s$^{-1}$]")
    tag = ("fixed divisor" if col == COL_FIXED else "OLD per-level divisor -- does NOT close") \
          + (", symlog" if scale != "lin" else ", linear")
    fig.suptitle(f"Meridional mass streamfunction  Psi(lat, z)  -  {OUT}  [{tag}]", y=0.98)
    sfx = ("" if col == COL_FIXED else "_olddiv") + ("" if scale != "lin" else "_lin")
    out = f"{OUT}/streamfn_lat_z_iter{iters[0]}_vs_{iters[1]}{sfx}.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)

    thr    = 0.25 * np.nanmax(np.abs(data[0][2][col]))
    for it, (lat, z, P) in zip(iters, data):
        strong = np.abs(P[col]) > thr
        zmax = z[np.where(strong.any(axis=1))[0].max()] if strong.any() else 0.0
        print(f"  iter {it:3d}: peak|Psi|={np.nanmax(np.abs(P[col])):6.1f}e9 kg/s   "
              f"overturning reaches up to {zmax/1000.0:4.2f} km (|Psi|>25% of iter-{iters[0]} peak)")
