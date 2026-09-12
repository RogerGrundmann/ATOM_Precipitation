#!/bin/bash
# _VW=0 AND THE RELEASED CEILING TOGETHER, 600 FROM SCRATCH.
#
# ONE VARIABLE AND ONE BINARY against output_pdc600: ATM_RADIAL_SHAPIRO_STRENGTH_VW 0.25 -> 0.
# Everything else identical, including ATM_PDYN_CEILING=50 and the DEFAULT ATM_BUOY_CONSISTENT=1.
#
# WHY THIS PAIRING AND NOT vw0. output_vw0 is also _VW=0 at 600 from scratch, but it carries
# ATM_BUOY_CONSISTENT=0, where max|p_dyn| is 4.3e-03 -- four orders below even the shipped 3.0
# ceiling -- so releasing the ceiling there is a guaranteed exact null and measures nothing.
# The ceiling only exists as a constraint on the =1 branch, which is the default.
#
# WHAT EACH HALF BOUGHT SEPARATELY (both measured 2026-09-12):
#   _VW=0            freezes the cells' AMPLITUDE (decay 20.3 % -> 1.7 % at 15N) and their FORM
#                    (normalised profile identical to three decimals over 380 iterations)
#   PDYN_CEILING=50  lets the pressure exceed 3.0 -- pdc600 reached 5.15 by iteration 560 with
#                    ZERO cells clipped where tr600d was truncating 4.8 % of the grid -- and
#                    raised the buoyancy cancellation 0.204 -> 0.291 on the 80-iteration arms
#
# ⚠ THIS IS THE LEAST-DAMPED CONFIGURATION THIS TREE HAS RUN: the v/w radial filter off, the
# pressure clamp 17x looser, and the buoyancy flip on (which drives the 51x radial runaway).
# u's OWN filter stays at full strength (ATM_RADIAL_SHAPIRO_STRENGTH=1.0 default), so the CFL
# guard the filter exists for is intact -- the hazard and the damage are in different components.
#
# EXPECTED, recorded before the run so it cannot be rationalised afterwards:
#   - the cells should keep their FORM (that is _VW=0's measured effect and the ceiling does not
#     enter the v equation), so the interesting question is NOT the cells but max|u|;
#   - max|u| is the open one. Two opposing arguments: the filter was removing radial-branch
#     momentum (so removing it could make the runaway worse), but _VW does not touch u's own
#     filter and the runaway is radial, so it may be untouched. pdc600 is the control.
#   - a NaN is a legitimate outcome and locates the load-bearing damping.
# WATCH: max|u|, max w_u against the +-100 clamp, cells at 20/200/400/600, CEILING BINDING,
# and iterations 155 / 357 / 483.
set -u
cd "$(dirname "$0")"
# start only when the pdc600 partner has finished, so the two do not share the 24 threads
while ! grep -q "exit" pdc600_driver.log 2>/dev/null; do sleep 30; done
export OMP_NUM_THREADS=24
echo "=== [vw0pdc] _VW=0 + PDYN_CEILING=50, nm=600 from scratch  $(date +%H:%M:%S) ==="
env ATM_CELLS_FROM_PSI=1 ATM_RADIAL_SHAPIRO_STRENGTH_VW=0 ATM_V_MASSBAL_STRIDE=1 \
    ATM_CELL_ROT_DIAG=1 ATM_TROPO_INDEX_FIX=1 ATM_UBUD_BALANCE=1 ATM_PDYN_CEILING=50.0 \
    ../cli/atm config_vw0pdc.xml > vw0pdc.log 2>&1
echo "  exit $?  nan $(grep -ci nan vw0pdc.log)  $(date +%H:%M:%S)"
