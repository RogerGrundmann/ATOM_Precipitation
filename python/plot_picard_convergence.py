#!/usr/bin/env python
# Convergence comparison for the atm<->ocean SST Picard loop:
# warm-start (ocean iterations accumulate) vs plain (from-scratch each round),
# for the hydrosphere (ocean) and the atmosphere.
#   * KE-drift %/window  = circulation convergence (converged when < 1%)
#   * mean T (K)         = thermal state
# Reads convergence.csv (atm) and convergence_hyd.csv (hyd) from each round dir.
import os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NM = 400
BLUE, VERM, INK, MUTED = "#0072B2", "#D55E00", "#1a1a1a", "#8a8a8a"  # CVD-safe (validated)
RUNS = {
    "warm-start": dict(base="output_picard_100Ma_ws",   rounds=4, color=BLUE),
    "plain":      dict(base="output_picard_100Ma",       rounds=3, color=VERM),
}
MODELS = {"hyd": ("convergence_hyd.csv", "Ocean (hydrosphere)"),
          "atm": ("convergence.csv",     "Atmosphere")}


def series(base, rounds, fname):
    """Return (x_cumulative_iter, drift_KE_pct, mean_T_K). Ocean warm-start rows
    already carry a cumulative iter (total_iter_count); everything that restarts
    each round is offset by round*NM so the x-axis is total effort spent."""
    # columns: iter, mean_T_K, drift_T_pct, drift_T_K, mean_KE_m2s2, drift_KE_pct, converged
    xs, ke, mT = [], [], []
    for r in range(rounds):
        p = os.path.join(base, "round_%d" % r, fname)
        if not os.path.exists(p):
            continue
        with open(p) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("iter"):     # skip blanks / any header (warm-start rounds have none)
                    continue
                c = line.split(",")
                it = int(c[0])
                x = it if it > NM else r * NM + it          # unified cumulative axis
                xs.append(x)
                mT.append(float(c[1]))
                ke.append(float(c[5]))
    return xs, ke, mT


fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.0), sharex=True)
fig.suptitle("atm↔ocean SST Picard loop — convergence: warm-start vs plain (Ma = 100)",
             fontsize=13, fontweight="bold", y=0.98)

for col, (mkey, (fname, mlabel)) in enumerate(MODELS.items()):
    a_ke, a_T = ax[0][col], ax[1][col]
    for rname, cfg in RUNS.items():
        xs, ke, mT = series(cfg["base"], cfg["rounds"], fname)
        # KE-drift: drop the window-warmup rows (drift==0 while the trailing window fills)
        xk = [x for x, v in zip(xs, ke) if v > 0]
        yk = [v for v in ke if v > 0]
        a_ke.plot(xk, yk, color=cfg["color"], lw=2, label=rname)
        a_T.plot(xs, mT, color=cfg["color"], lw=2, label=rname)
        # direct label at each curve's end (secondary encoding, not color-alone)
        if xk:
            a_ke.annotate(rname, (xk[-1], yk[-1]), color=cfg["color"], fontsize=9,
                          fontweight="bold", xytext=(5, 0), textcoords="offset points", va="center")
    # 1% converged threshold on the KE-drift panels
    a_ke.axhline(1.0, color=MUTED, lw=1, ls="--")
    a_ke.annotate("1% converged", (a_ke.get_xlim()[0], 1.0), color=MUTED, fontsize=8,
                  xytext=(4, 3), textcoords="offset points")
    a_ke.set_title(mlabel, fontsize=11, color=INK)
    a_ke.set_ylim(bottom=0)
    a_ke.set_ylabel("KE drift  [% / window]", fontsize=9)
    a_T.set_ylabel("mean T  [K]", fontsize=9)
    a_T.set_xlabel("cumulative iteration (effort)", fontsize=9)
    a_T.ticklabel_format(axis="y", useOffset=False)         # no +2.595e2 offset box
    if mkey == "atm":
        a_ke.annotate("(both runs overlap — atm always from scratch)", (0.5, 0.9),
                      xycoords="axes fraction", ha="center", fontsize=8, color=MUTED)
    for axi in (a_ke, a_T):
        for s in ("top", "right"):
            axi.spines[s].set_visible(False)
        axi.grid(True, color="#e6e6e6", lw=0.6)
        axi.set_axisbelow(True)
        axi.tick_params(labelsize=8)
    # round-boundary ticks
    for r in range(1, max(c["rounds"] for c in RUNS.values())):
        a_ke.axvline(r * NM, color="#f0f0f0", lw=1, zorder=0)
        a_T.axvline(r * NM, color="#f0f0f0", lw=1, zorder=0)

ax[0][0].legend(frameon=False, fontsize=9, loc="upper right")
fig.text(0.5, 0.012,
         "Ocean KE-drift falls 5.5→1.8% as iterations accumulate (warm-start) but stays stuck "
         "at 5.5% every round when the ocean restarts (plain).",
         ha="center", fontsize=9, color=MUTED)
fig.tight_layout(rect=[0, 0.035, 1, 0.96])
out = "picard_convergence.png"
fig.savefig(out, dpi=140)
print("wrote", out)
