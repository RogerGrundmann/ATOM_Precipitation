#!/bin/bash
# ATM_BUOY_TREF, THE MIDDLE ARM, RUN AT LAST (2026-09-14).
#
# (b) alone: divide the Boussinesq anomaly by t_ref_level[i] -- the PHYSICAL reference
# temperature -- instead of leaving it divided by t_0 = 273.15 K, the NON-DIMENSIONALISATION
# constant. Unreachable from 2026-09-09 to 2026-09-14 because ATM_BUOY_CONSISTENT implied it and
# was the default; reachable again since the revert (commit 0552b9a).
#
# 600 -> 700 from output_vw0/atm_restart_0Ma_600.bin, which is a 600-from-scratch run on EXACTLY
# the post-revert default configuration (its banner: SATADJ_PHASE=1 _VW=0 CELLS_FROM_PSI=1
# V_MASSBAL_STRIDE=1 METRIC_SIN_FLOOR=0.26 TROPO_INDEX_FIX=1 BUOY_CONSISTENT=0 PDYN_CEILING=3.0).
# Both arms seeded from the identical .bin, md5 b0b7f85f.
#
# WHY A RESTART AND WHY 100 ITERATIONS. buoyancy_ramp = min(1, total_iter/300), so from 600 the
# term is at FULL strength for the whole window -- which the README's own note says a 4-iteration
# A/B could not do (ramp 0.013 there). And ubud_* are filled only when do_vbudget fires, and the
# FIRST checkpoint after a restart yields zeros, so an arm shorter than ~40 iterations would print
# a full table of zeros. 100 gives four usable budget prints.
#
# PREDICTION, RECORDED BEFORE THE RUN.
#  1. The multiplier is 1/t_ref_level[i], computed from vw0's own iteration-600 field:
#     0.955 at the ground, 0.969 at 1.4 km, 1.005 at 3.3 km, 1.097 at 7.4 km, 1.260 at the lid.
#     So this is NOT the "~5 % correction" CLAUDE.md and the README call it -- that is the SURFACE
#     value. It is a 32 % height-dependent distortion THAT CHANGES SIGN AT 3.3 km: weaker below,
#     up to +26 % stronger above.
#  2. Therefore rms ubud_buoy should move by the amplitude-weighted mean multiplier. The
#     temperature anomaly is largest low down, so I expect it DOWN by a few per cent, not up.
#  3. Everything else NULL. Shipped ubud_buoy is 2.46e-05 against ubud_pgf 8.23e-02 = 0.03 %, so
#     a +-5..26 % change to it moves rhs_u by at most 0.008 %. max|u|, Psi, precipitation, r and
#     the bands should be unchanged to the printed digits.
#  If 3 fails -- if anything moves -- then the buoyancy is NOT inert on this branch and the
#  0.03 % figure, which was measured at iteration 700 of a different configuration, is wrong.
#
# PARALLEL, 12 THREADS EACH (see verify_lib.sh): two independent processes, and both arms use the
# SAME thread count, which is what comparability requires -- not that the count be 24.
set -u
cd "$(dirname "$0")"
export OMP_NUM_THREADS=12
rm -f BTREF_DONE
echo "=== btref 600->700, both arms in parallel, 12 threads each  $(date +%H:%M:%S)"
( env ATM_UBUD_BALANCE=1                   ../cli/atm config_bt_ctl.xml > bt_ctl.log 2>&1
  echo "  bt_ctl exit $?  nan $(grep -ci 'nan' bt_ctl.log)  $(date +%H:%M:%S)" ) &
( env ATM_UBUD_BALANCE=1 ATM_BUOY_TREF=1   ../cli/atm config_bt_on.xml  > bt_on.log  2>&1
  echo "  bt_on  exit $?  nan $(grep -ci 'nan' bt_on.log)   $(date +%H:%M:%S)" ) &
wait
echo "=== btref done $(date +%H:%M:%S)"
touch BTREF_DONE
