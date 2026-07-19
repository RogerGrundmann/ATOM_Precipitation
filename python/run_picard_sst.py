#!/usr/bin/env python
# ---------------------------------------------------------------------------
# atm<->hyd SST Picard outer loop (scaffold).
#
# Fixed-point iteration on the ocean sea-surface temperature. The frozen Scotese
# zonal parabola is the round-0 initial guess; each round re-runs BOTH models and
# blends the hydrosphere's SST back into the atmosphere until the SST stops moving.
#
#   round 0:  ATM (alpha=0, no SST file yet) -> winds
#             HYD (reads winds)              -> SST_0
#   round r:  copy SST_{r-1} into round dir
#             ATM (sst_coupling_alpha=alpha, reads SST_{r-1}, blends into t.x[0]) -> winds
#             HYD (reads winds)                                                   -> SST_r
#             convergence:  max|SST_r - SST_{r-1}|  (over ocean, reported in K)
#
# The reverse channel is the C++ pair HydrosphereSSTTransfer (writes
# <stem>_Transfer_Hyd_SST_<iter>.vwtp) / read_Hydrosphere_SST (blends it in). Each
# round runs from scratch, so total_iter_count restarts and the SST stamp repeats;
# rounds are therefore isolated in their own sub-directories and the previous
# round's SST is copied forward explicitly (below), rather than relying on the
# in-dir "latest snapshot" auto-select.
#
# Safety: round 0 is EXACTLY the existing one-way chain (alpha forced to 0, no SST
# file present). The C++ blend is ocean-only, finite-checked, SST-clamped and
# ocean-mean-anchored, so it moves SST structure without shifting ocean energy.
#
# Usage:
#   python run_picard_sst.py <Ma> [--rounds N] [--alpha A] [--nm ITERS]
#                                  [--tol K] [--base DIR]
# Example:
#   python run_picard_sst.py 100 --rounds 4 --alpha 0.4 --nm 400
# ---------------------------------------------------------------------------
import sys, os, glob, shutil, argparse
from pyatom import Atmosphere, Hydrosphere

T_0 = 273.15  # non-dim temperature reference: T_C = (t_nd - 1) * T_0, so dT_C = dt_nd * T_0


def parse_args():
    p = argparse.ArgumentParser(description="atm<->hyd SST Picard outer loop")
    p.add_argument("Ma", type=int, help="paleo time slice in Ma")
    p.add_argument("--rounds", type=int, default=4, help="number of Picard rounds (incl. round 0)")
    p.add_argument("--alpha", type=float, default=0.4, help="atm-side SST blend fraction for rounds >= 1")
    p.add_argument("--hyd-relax", type=float, default=0.5,
                   help="hydrosphere sst_relax_alpha (surface flux-BC strength). MUST be < 1 for a "
                        "non-trivial loop: at 1.0 the ocean surface is hard re-pinned to the atm SST, "
                        "so the reverse channel just returns the same field and the loop converges "
                        "trivially. 0 = unforced cold-collapse (do not use).")
    p.add_argument("--nm", type=int, default=400, help="iterations per atm/hyd run")
    p.add_argument("--checkpoint", type=int, default=100,
                   help="VTK/panorama printout stride within each run (also the report-only "
                        "convergence-monitor sampling cadence; convergence.csv is written per round dir)")
    p.add_argument("--tol", type=float, default=0.05, help="Picard convergence tol on max|dSST| in K")
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


def max_dsst_kelvin(path_a, path_b):
    """max|a - b| over cells, in K. Land is 0 in both files so contributes 0."""
    a, b = read_sst(path_a), read_sst(path_b)
    n = min(len(a), len(b))
    return max((abs(a[i] - b[i]) for i in range(n)), default=0.0) * T_0


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


def run_hyd(Ma, outdir, nm, hyd_relax, ckpt):
    h = Hydrosphere()
    h.output_path = outdir.encode()
    h.input_path = outdir.encode()
    h.nm = nm
    h.checkpoint = ckpt                        # VTK/MinMax printouts every ckpt iters
    h.panorama_print = ckpt                    # panorama every ckpt iters
    h.sst_relax_alpha = hyd_relax             # < 1 => Haney flux BC, ocean carries its own SST anomaly
    print("    [HYD] Ma=%d nm=%d ckpt=%d sst_relax_alpha=%.3f -> %s" % (Ma, nm, ckpt, hyd_relax, outdir), flush=True)
    h.run_time_slice(Ma)


def main():
    args = parse_args()
    base = args.base or ("output_picard_%dMa" % args.Ma)
    os.makedirs(base, exist_ok=True)
    print("### SST Picard loop: Ma=%d rounds=%d alpha=%.3f nm=%d tol=%.3fK base=%s ###"
          % (args.Ma, args.rounds, args.alpha, args.nm, args.tol, base), flush=True)

    prev_sst = None
    for r in range(args.rounds):
        rdir = os.path.join(base, "round_%d" % r) + "/"
        os.makedirs(rdir, exist_ok=True)
        alpha = 0.0 if r == 0 else args.alpha

        # Carry the previous round's SST into this round's dir so the atm reader
        # (which scans its own output_path) can find and blend it.
        if r >= 1 and prev_sst is not None:
            shutil.copy(prev_sst, rdir)

        print("=== round %d/%d (alpha=%.3f) ===" % (r, args.rounds - 1, alpha), flush=True)
        run_atm(args.Ma, rdir, args.nm, alpha, args.checkpoint)
        run_hyd(args.Ma, rdir, args.nm, args.hyd_relax, args.checkpoint)

        cur_sst = latest_sst_file(rdir)
        if cur_sst is None:
            print("    !! no SST snapshot produced in %s — hyd run failed? stopping." % rdir, flush=True)
            return 2

        if prev_sst is not None:
            d = max_dsst_kelvin(cur_sst, prev_sst)
            print("    >> round %d: max|dSST| = %.4f K (tol %.4f K)" % (r, d, args.tol), flush=True)
            if d < args.tol:
                print("### CONVERGED at round %d (max|dSST|=%.4f K < %.4f K) ###" % (r, d, args.tol), flush=True)
                return 0
        prev_sst = cur_sst

    print("### reached round limit (%d) without hitting tol ###" % (args.rounds - 1), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
