#!/usr/bin/env python3
"""Reconstruct the hydro convergence-monitor metrics (ocean-mean T, mean KE)
directly from the saved restart .bin checkpoints — no model re-run required.
Prints a per-checkpoint table and writes a two-panel convergence plot.

Why from bins: the hydro convergence monitor (cHydrosphereModel.cpp, writes
convergence_hyd.csv) was added AFTER the Ma=100/200 paleo 1000-iter runs, so
those runs never produced the CSV. Every checkpoint holds the full state, so
the same two integral metrics can be computed offline.

Bin layout (FileIO_Hyd.cpp save_state): int32 hdr[5]={magic 0x4F434D31,im,jm,km,
total_iter_count} then 17 float64 arrays in restart_arrays() order
[t,u,v,w,c,tn,un,vn,wn,cn,p_dyn,p_hydro,tke,dis,tken,disn,nue], each im*jm*km,
layout [i][j][k]. t=idx0, u=1, v=2, w=3.

Land mask (h) is NOT in the restart, so land is detected as cells whose
velocity is exactly 0 in EVERY checkpoint (bcSolidGround zeroes land velocity
each iter). Ocean-mean is sin(colat) area-weighted; the drift trend is
weighting-independent.

Usage: python3 analyze_hyd_convergence.py   (run from the python/ dir)
"""
import numpy as np, glob, re, os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

t_0, u_0 = 273.15, 0.24                       # cHydrosphereDefaults.cpp.inc
IM, JM, KM = 41, 181, 361
NCELL = IM * JM * KM
HDR = 5 * 4
BLUE, ORANGE = "#2a78d6", "#eb6834"          # CVD-validated categorical pair
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e6e6e2", "#fcfcfb"

RUNS = [(100, "output_100Ma_paleo"), (200, "output_200Ma_paleo")]


def load(fn, which):
    return np.fromfile(fn, dtype="<f8", count=NCELL,
                       offset=HDR + which * NCELL * 8).reshape(IM, JM, KM)


def series(ddir, Ma):
    bins = sorted(glob.glob(f"{ddir}/hyd_restart_{Ma}Ma_*.bin"),
                  key=lambda f: int(re.search(r"_(\d+)\.bin$", f).group(1)))
    if not bins:
        return None
    iters = [int(re.search(r"_(\d+)\.bin$", f).group(1)) for f in bins]
    land = np.ones((IM, JM, KM), bool)
    for f in bins:
        land &= (load(f, 1) == 0) & (load(f, 2) == 0) & (load(f, 3) == 0)
    water = ~land
    wt = np.broadcast_to(np.sin(np.linspace(0, np.pi, JM))[None, :, None],
                         (IM, JM, KM)) * water
    ws = wt.sum()
    T, KE = [], []
    for f in bins:
        T.append((load(f, 0) * wt).sum() / ws * t_0)
        KE.append((0.5 * (load(f, 1)**2 + load(f, 2)**2 + load(f, 3)**2)
                   * wt).sum() / ws * u_0 * u_0)
    return np.array(iters), np.array(T), np.array(KE), water.mean()


def main():
    data = {}
    for Ma, ddir in RUNS:
        s = series(ddir, Ma)
        if s is None:
            print(f"(no bins in {ddir})"); continue
        it, T, KE, frac = s
        data[Ma] = s
        print(f"\n=== Ma={Ma}  ({ddir})  ocean frac={frac:.3f} ===")
        print(f"{'iter':>6} {'meanT_K':>10} {'dT_K':>9} {'meanKE_m2s2':>13} {'dKE%':>8}")
        pT = pKE = None
        for i, T_i, KE_i in zip(it, T, KE):
            dT = f"{T_i-pT:+.4f}" if pT is not None else "   --"
            dKE = f"{100*(KE_i-pKE)/max(1e-12,KE_i):+.2f}" if pKE is not None else "  --"
            print(f"{i:>6} {T_i:>10.4f} {dT:>9} {KE_i:>13.6e} {dKE:>8}")
            pT, pKE = T_i, KE_i

    if not data:
        return
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK})
    fig, (axT, axK) = plt.subplots(1, 2, figsize=(11, 4.4), facecolor=SURF)
    for ax in (axT, axK):
        ax.set_facecolor(SURF); ax.grid(True, color=GRID, lw=.8, zorder=0)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        ax.set_xlabel("iteration")
    color = {100: BLUE, 200: ORANGE}
    for Ma in sorted(data):
        it, T, KE, _ = data[Ma]
        axT.plot(it, T, color=color[Ma], lw=2, marker="o", ms=5, zorder=3)
        axK.plot(it, KE * 1e4, color=color[Ma], lw=2, marker="o", ms=5, zorder=3)
        for ax, y in ((axT, T[-1]), (axK, KE[-1] * 1e4)):
            ax.annotate(f"Ma={Ma}", (it[-1], y), color=color[Ma], fontweight="bold",
                        xytext=(8, 0), textcoords="offset points", va="center", fontsize=10)
    axT.set_title("Ocean-mean temperature", fontsize=12, fontweight="bold", loc="left", pad=10)
    axT.set_ylabel("mean T  [K]")
    axK.set_title("Ocean-mean kinetic energy", fontsize=12, fontweight="bold", loc="left", pad=10)
    axK.set_ylabel("mean KE  [10⁻⁴ m²/s²]")
    fig.suptitle("Hydrosphere convergence from restart checkpoints",
                 fontsize=13, fontweight="bold", x=.02, ha="left", y=.99)
    fig.tight_layout(rect=[0, 0, 1, .95])
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "hyd_convergence_1000.png")
    fig.savefig(out, dpi=150, facecolor=SURF)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
