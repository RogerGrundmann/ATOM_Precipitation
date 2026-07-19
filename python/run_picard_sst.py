#!/usr/bin/env python
# ---------------------------------------------------------------------------
# atm<->hyd SST Picard outer loop, with optional Anderson acceleration.
#
# Fixed-point iteration on the ocean sea-surface temperature. The frozen Scotese
# zonal parabola is the round-0 initial guess; each round re-runs BOTH models and
# blends the hydrosphere's SST back into the atmosphere until the SST stops moving.
#
#   round 0:  ATM (alpha=0, no SST file yet) -> winds ;  HYD (reads winds) -> SST_0 = g_0
#   round r:  form the input field u_r (below) and write it where the atm reads it
#             ATM (sst_coupling_alpha=alpha, reads u_r, blends into t.x[0]) -> winds
#             HYD (reads winds)                                             -> g_r
#             fixed-point residual:  max|g_r - u_r|  (over ocean, reported in K)
#
# Define the round map G by g_r = G(u_r): "atm blends u_r, runs; hyd runs; returns SST".
# The fixed point u* = g* (the round returns the SST it was given) is a self-consistent
# atm<->ocean state. Two ways to pick u_r:
#
#   * PLAIN PICARD (--anderson 0, default): u_r = g_{r-1}. This is exactly the loop
#     validated at Ma=100 (round-delta contraction ~9x/round). numpy is not used.
#
#   * ANDERSON(m) (--anderson m>0): u_{k+1} = u_k + beta*f_k - (dU + beta*dF) gamma,
#     f_k = g_k - u_k, gamma = argmin || f_k - dF gamma ||_2, where dF/dU are the last
#     m residual/input differences (Walker & Ni 2011). It mixes a short history of past
#     iterates to cancel the leading error modes, so a well-behaved loop reaches the
#     fixed point in fewer rounds. With m=0 or too little history it falls back to the
#     damped Picard step u_k + beta*f_k (== plain Picard at beta=1). Anderson genuinely
#     engages once >=2 residuals exist (round 3+), so run --rounds 4-5 to benefit.
#     Overshoot is safe: the C++ read_Hydrosphere_SST clamps SST to [-1.8,40] C and
#     re-anchors the ocean-mean before blending.
#
# The reverse channel is the C++ pair HydrosphereSSTTransfer (writes
# <stem>_Transfer_Hyd_SST_<iter>.vwtp) / read_Hydrosphere_SST (blends it in). Each round
# runs from scratch (total_iter_count restarts, so the SST stamp repeats); rounds are
# isolated in their own sub-directories and the input field is placed there explicitly.
#
# WARM-START (--warm-start): by default every round re-spins the OCEAN from its IC, so
# only ~nm iters of ocean adjustment ever accumulate — and since the ocean is far from
# equilibrium at nm~400 (KE still drifting), the rounds are near-identical partial
# spin-downs and the loop converges to the fixed point of the *nm-iter* map, not the true
# ocean equilibrium. --warm-start instead resumes each round's hydrosphere from the
# previous round's restart .bin, so ocean iterations ACCUMULATE (R rounds -> R*nm effective
# ocean iters, progressively approaching equilibrium) -- the asynchronous Manabe-Bryan
# style of paleo coupling. This is safe because restart_arrays() restores only the ocean
# interior (t,u,v,w,c,...); the wind forcing (v_wind/w_wind) and SST target (t_surf_fix)
# are re-captured at init from THIS round's atmosphere, so the ocean carries over while the
# forcing is refreshed. The ATMOSPHERE still runs from scratch each round (it is the fast
# component and must re-init to blend the new SST). Needs nm a multiple of --checkpoint so
# the final iter is checkpointed.
#
# Safety: round 0 is EXACTLY the existing one-way chain (alpha forced to 0, no SST file).
#
# Usage:
#   python run_picard_sst.py <Ma> [--rounds N] [--alpha A] [--hyd-relax R] [--nm ITERS]
#                                  [--anderson M] [--anderson-beta B]
#                                  [--checkpoint C] [--tol K] [--base DIR]
# Example:
#   python run_picard_sst.py 100 --rounds 5 --alpha 0.4 --hyd-relax 0.5 --anderson 2
# ---------------------------------------------------------------------------
import sys, os, glob, shutil, argparse
from pyatom import Atmosphere, Hydrosphere

T_0 = 273.15  # non-dim temperature reference: T_C = (t_nd - 1) * T_0, so dT_C = dt_nd * T_0


def parse_args():
    p = argparse.ArgumentParser(description="atm<->hyd SST Picard outer loop (+ Anderson accel)")
    p.add_argument("Ma", type=int, help="paleo time slice in Ma")
    p.add_argument("--rounds", type=int, default=4, help="number of Picard rounds (incl. round 0)")
    p.add_argument("--alpha", type=float, default=0.4, help="atm-side SST blend fraction for rounds >= 1")
    p.add_argument("--hyd-relax", type=float, default=0.5,
                   help="hydrosphere sst_relax_alpha (surface flux-BC strength). MUST be < 1 for a "
                        "non-trivial loop: at 1.0 the ocean surface is hard re-pinned to the atm SST, "
                        "so the reverse channel just returns the same field and the loop converges "
                        "trivially. 0 = unforced cold-collapse (do not use).")
    p.add_argument("--anderson", type=int, default=0,
                   help="Anderson acceleration history depth m (0 = plain Picard). m~2-3 is typical; "
                        "engages from round 3 on. Needs numpy.")
    p.add_argument("--anderson-beta", type=float, default=1.0,
                   help="Anderson mixing/damping (1.0 = none; <1 adds under-relaxation for robustness).")
    p.add_argument("--warm-start", action="store_true",
                   help="resume each round's hydrosphere from the previous round's restart .bin, so "
                        "ocean iterations accumulate across rounds (R*nm effective) instead of "
                        "re-spinning from the IC every round. The atmosphere still runs from scratch.")
    p.add_argument("--nm", type=int, default=400, help="iterations per atm/hyd run")
    p.add_argument("--checkpoint", type=int, default=100,
                   help="VTK/panorama printout stride within each run (also the report-only "
                        "convergence-monitor sampling cadence; convergence.csv is written per round dir)")
    p.add_argument("--tol", type=float, default=0.05,
                   help="stop when the fixed-point residual max|g_r - u_r| drops below this (K)")
    p.add_argument("--base", default=None, help="base output dir (default output_picard_<Ma>Ma)")
    return p.parse_args()


def latest_sst_file(d):
    """Newest hydrosphere SST snapshot in dir d, or None."""
    files = glob.glob(os.path.join(d, "*_Transfer_Hyd_SST_*.vwtp"))
    if not files:
        return None
    def iter_of(f):
        stem = f[:-len(".vwtp")]
        return int(stem.rsplit("_", 1)[1])
    return max(files, key=iter_of)


def latest_restart(d, Ma):
    """Newest hydrosphere restart checkpoint in dir d as (path, iter), or (None, -1)."""
    files = glob.glob(os.path.join(d, "hyd_restart_%dMa_*.bin" % Ma))
    if not files:
        return None, -1
    def iter_of(f):
        return int(os.path.basename(f)[:-len(".bin")].rsplit("_", 1)[1])
    f = max(files, key=iter_of)
    return f, iter_of(f)


def read_sst(path):
    """Return a flat list of non-dim SST values (skips the '# iter_n' headline)."""
    vals = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            vals.append(float(line.split()[0]))
    return vals


def write_sst(path, vals):
    """Write a non-dim SST field in the format read_Hydrosphere_SST expects
    (one value per row, a leading '# iter_n' headline). Matches the C++ writer's
    6-decimal format so the atm reads a driver-supplied field just like a hyd one."""
    with open(path, "w") as f:
        f.write("# iter_n = picard_input\n")
        f.write("\n".join("%.6f" % v for v in vals))
        f.write("\n")


def max_resid_kelvin(g_vals, u_vals):
    """Fixed-point residual max|g - u| in K. Land is 0 in both, so it drops out."""
    n = min(len(g_vals), len(u_vals))
    return max((abs(g_vals[i] - u_vals[i]) for i in range(n)), default=0.0) * T_0


def anderson_next(U, G, m, beta):
    """Anderson(m) next input field. U, G are lists of value-lists (inputs actually
    used, and the hyd outputs they produced), most recent last. Returns a value-list.
    Walker & Ni (2011): with f_i = g_i - u_i,
        gamma = argmin || f_k - dF gamma ||,   dF = [f_{i+1}-f_i] over the last min(m,k) steps
        u_next = u_k + beta*f_k - (dU + beta*dF) gamma
    Falls back to the damped Picard step u_k + beta*f_k when there is no history."""
    import numpy as np
    Ua = [np.asarray(u) for u in U]
    Ga = [np.asarray(g) for g in G]
    F = [g - u for u, g in zip(Ua, Ga)]
    k = len(F) - 1
    u_k, f_k = Ua[k], F[k]
    hist = min(m, k)
    if hist == 0:
        return (u_k + beta * f_k).tolist()
    dF = np.column_stack([F[i + 1] - F[i] for i in range(k - hist, k)])
    dU = np.column_stack([Ua[i + 1] - Ua[i] for i in range(k - hist, k)])
    gamma, *_ = np.linalg.lstsq(dF, f_k, rcond=None)
    return (u_k + beta * f_k - (dU + beta * dF) @ gamma).tolist()


def run_atm(Ma, outdir, nm, alpha, ckpt):
    a = Atmosphere()
    a.output_path = outdir.encode()          # std::string binding wants bytes in py3
    a.nm = nm
    a.checkpoint = ckpt                       # VTK/MinMax printouts every ckpt iters
    a.panorama_print = ckpt                   # panorama every ckpt iters
    a.sst_coupling_alpha = alpha              # 0.0 in round 0 -> no SST read
    a.hyd_sst_iter = -1                       # latest snapshot in outdir
    print("    [ATM] Ma=%d nm=%d ckpt=%d alpha=%.3f -> %s" % (Ma, nm, ckpt, alpha, outdir), flush=True)
    a.run_time_slice(Ma)


def run_hyd(Ma, outdir, nm, hyd_relax, ckpt, restart_from_iter=-1):
    h = Hydrosphere()
    h.output_path = outdir.encode()
    h.input_path = outdir.encode()
    h.nm = nm                                  # ADDITIONAL iters on a restart (total_iter_count accumulates)
    h.checkpoint = ckpt                        # VTK/MinMax printouts every ckpt iters
    h.panorama_print = ckpt                    # panorama + restart .bin every ckpt iters
    h.sst_relax_alpha = hyd_relax             # < 1 => Haney flux BC, ocean carries its own SST anomaly
    h.restart_from_iter = restart_from_iter    # >=0 => warm-start the ocean from that checkpoint
    tag = (" restart@%d" % restart_from_iter) if restart_from_iter >= 0 else ""
    print("    [HYD] Ma=%d nm=%d ckpt=%d sst_relax_alpha=%.3f%s -> %s"
          % (Ma, nm, ckpt, hyd_relax, tag, outdir), flush=True)
    h.run_time_slice(Ma)


def main():
    args = parse_args()
    base = args.base or ("output_picard_%dMa" % args.Ma)
    os.makedirs(base, exist_ok=True)
    accel = "Anderson(%d, beta=%.2f)" % (args.anderson, args.anderson_beta) if args.anderson > 0 else "plain Picard"
    ocean = "warm-start" if args.warm_start else "from-scratch/round"
    print("### SST Picard loop: Ma=%d rounds=%d alpha=%.3f hyd_relax=%.3f nm=%d tol=%.3fK %s ocean=%s base=%s ###"
          % (args.Ma, args.rounds, args.alpha, args.hyd_relax, args.nm, args.tol, accel, ocean, base), flush=True)

    prev_sst = None       # path to g_{r-1}'s SST file
    prev_rdir = None      # previous round's dir (source of the ocean restart .bin)
    U_hist = []           # value-lists: input fields actually used (rounds >= 1)
    G_hist = []           # value-lists: hyd outputs g_r (aligned with U_hist)
    for r in range(args.rounds):
        rdir = os.path.join(base, "round_%d" % r) + "/"
        os.makedirs(rdir, exist_ok=True)
        alpha = 0.0 if r == 0 else args.alpha

        # Build the input field u_r and place it where round r's atm will read it.
        if r >= 1:
            name = os.path.basename(prev_sst)   # keep the stem/iter so the atm's scan selects it
            if args.anderson > 0 and len(U_hist) >= 1:
                u_vals = anderson_next(U_hist, G_hist, args.anderson, args.anderson_beta)
                write_sst(os.path.join(rdir, name), u_vals)
            else:
                shutil.copy(prev_sst, rdir)     # plain Picard: input = previous output (byte-identical)
                u_vals = read_sst(prev_sst)
            U_hist.append(u_vals)

        # Warm-start: carry the previous round's ocean restart into this round's dir so the
        # hydrosphere resumes from it (ocean iters accumulate). The atmosphere stays fresh.
        restart_iter = -1
        if args.warm_start and r >= 1 and prev_rdir is not None:
            src_bin, src_iter = latest_restart(prev_rdir, args.Ma)
            if src_bin is not None:
                shutil.copy(src_bin, rdir)
                restart_iter = src_iter
                print("    [warm-start] ocean resumes from %s (total_iter %d)"
                      % (os.path.basename(src_bin), src_iter), flush=True)
            else:
                print("    [warm-start] no restart .bin in %s — ocean starts fresh this round" % prev_rdir, flush=True)

        print("=== round %d/%d (alpha=%.3f) ===" % (r, args.rounds - 1, alpha), flush=True)
        run_atm(args.Ma, rdir, args.nm, alpha, args.checkpoint)
        run_hyd(args.Ma, rdir, args.nm, args.hyd_relax, args.checkpoint, restart_from_iter=restart_iter)

        cur_sst = latest_sst_file(rdir)
        if cur_sst is None:
            print("    !! no SST snapshot produced in %s — hyd run failed? stopping." % rdir, flush=True)
            return 2
        g_vals = read_sst(cur_sst)

        if r >= 1:
            G_hist.append(g_vals)
            resid = max_resid_kelvin(g_vals, U_hist[-1])
            print("    >> round %d: max|dSST| (residual) = %.4f K (tol %.4f K)" % (r, resid, args.tol), flush=True)
            if resid < args.tol:
                print("### CONVERGED at round %d (residual=%.4f K < %.4f K) ###" % (r, resid, args.tol), flush=True)
                return 0
        prev_sst = cur_sst
        prev_rdir = rdir

    print("### reached round limit (%d) without hitting tol ###" % (args.rounds - 1), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
