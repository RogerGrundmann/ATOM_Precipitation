#!/bin/bash
# IS THE PRESSURE AMPLITUDE LIMITED BY THE p_dyn CEILING RATHER THAN BY THE SOLVER?
#
# Found in output_lsbal / output_ubal2, both already on disk and both printing it all along:
#   ubal2 (pointwise)  max|p_dyn| pre-clip 3.076  ->  5.24-5.35 % of ALL cells CLIPPED at 3.0
#   lsbal (line solve) max|p_dyn| pre-clip 3.415  ->  7.21-7.60 % CLIPPED
# and in tr600d (600 from scratch, current defaults) max|p_dyn| grows monotonically
# 7.9e-04 (iter 20) -> 2.971 and first binds at iteration ~400, reaching 4.8 % by 600.
# With ATM_BUOY_CONSISTENT=0 it never binds at all (vw0 4.3e-03, bcs0 4.2e-02, sp600ctl 1.7e-02).
#
# So CLAUDE.md's "neither clamp binds today, max|p_dyn| ~ 0.017, 176x below the ceiling" is
# a statement about the PRE-FLIP branch and stopped being true when ATM_BUOY_CONSISTENT
# became the default on 2026-09-09.
#
# CONSEQUENCE FOR YESTERDAY'S CONCLUSION: the line solve raises the pre-clip pressure 3.076
# -> 3.415 (it IS converging the column mode the body force excites) and the ceiling then
# truncates 42 % MORE cells, so the DELIVERED gradient goes DOWN, rms pgf 1.169 -> 1.116.
# "The amplitude barely moves, so the wide-gradient correction is the binding constraint"
# was therefore measured through a clamp and is not established.
#
# ONE VARIABLE PER ARM: ATM_PDYN_CEILING 3.0 -> 50.0, everything else as the partner run.
#   pdc_ls vs output_lsbal  (line solve)
#   pdc_pw vs output_ubal2  (pointwise; cli/atm with the knob unset is byte-identical to
#                            atm_rot600, established by the lb_a/lb_b check)
# 50 is above the ~40 non-dim a geostrophically balanced mid-latitude field needs, and 15x
# the shipped phase-dependent value.
#
# EXPECTED, recorded before the run so it cannot be rationalised afterwards:
#   - if rms pgf and the cancellation RISE, the clamp was the binding constraint and the
#     face-consistent correction is NOT the next thing to build;
#   - if they do not move, the amplitude limit is real and yesterday's reading stands;
#   - a NaN is a legitimate outcome: the ceiling was lowered to 3.0 for the iter-483 runaway,
#     so raising it may re-trigger it. That would itself be the answer.
# WATCH max|u| (the 51x runaway) and max w_u against the +-100 clamp.
set -u
cd "$(dirname "$0")"
export OMP_NUM_THREADS=24
B="ATM_CELLS_FROM_PSI=1 ATM_RADIAL_SHAPIRO_STRENGTH_VW=0.25 ATM_V_MASSBAL_STRIDE=1 ATM_TROPO_INDEX_FIX=1 ATM_UBUD_BALANCE=1 ATM_PDYN_CEILING=50.0"
echo "=== [pdc_ls] line solve + PDYN_CEILING=50, 600->680  $(date +%H:%M:%S) ==="
env $B ATM_PRESS_LINE_SOLVE=1 ../cli/atm config_pdc_ls.xml > pdc_ls.log 2>&1
echo "  exit $?  nan $(grep -ci nan pdc_ls.log)  $(date +%H:%M:%S)"
echo "=== [pdc_pw] pointwise  + PDYN_CEILING=50, 600->680  $(date +%H:%M:%S) ==="
env $B ../cli/atm config_pdc_pw.xml > pdc_pw.log 2>&1
echo "  exit $?  nan $(grep -ci nan pdc_pw.log)  $(date +%H:%M:%S)"
echo "=== queue done $(date +%H:%M:%S) ==="
