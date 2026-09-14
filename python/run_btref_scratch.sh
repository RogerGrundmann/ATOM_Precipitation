#!/bin/bash
# ATM_BUOY_TREF FROM SCRATCH TO 100, at the user's instruction: watch whether `u` is still
# growing and disturbing the cells (2026-09-14). Companion to run_btref.sh, which is the
# FULL-STRENGTH restart pair; this one is the TRAJECTORY pair.
#
# checkpoint = 10, so the streamfunction CSV and both momentum budgets land every 10 iterations
# and `max u-component` prints every iteration -- the growth curve, not two endpoints.
# ATM_CELL_ROT_DIAG=1 so the cells' ascent/descent boundary is reported as it evolves; it is the
# diagnostic that went silent past iteration ~60 under the runaway, printing "zonal-mean u does
# not change sign between 60 and 90 deg", and it is the direct read on `u` disturbing the cells.
#
# ⚠ THE RAMP IS THE LIMITATION OF THIS PAIR AND IT IS NOT SMALL. buoyancy_ramp = min(1, iter/300),
# so over 0..100 it runs 0 -> 0.33 and AVERAGES ~0.17: the term under test acts at about a sixth
# of full strength, which is exactly the under-powering the README's own note gives as the reason
# the first 4-iteration A/B of this knob measured nothing. So a NULL here is weak evidence and a
# non-null is strong evidence. The restart pair (600->700, ramp = 1.0 throughout) carries the
# full-strength measurement; read them together.
#
# EXPECTATION. On the post-revert default `u` DECAYS from scratch -- 0.0075 rms at iteration 20 to
# 0.0039 at 200 on the 87E slice -- so the control should NOT grow. The question this pair answers
# is whether ATM_BUOY_TREF re-introduces growth, which is plausible in DIRECTION because the knob
# multiplies the buoyancy by 1/t_ref_level[i]: 0.955 at the ground but 1.097 at 7.4 km and 1.260
# at the lid, i.e. it STRENGTHENS the force in exactly the upper half of the column where the
# BUOY_CONSISTENT runaway lived. If max|u| grows in bs_on and not in bs_ctl, that is the same
# mechanism at 1/6 strength and the knob should stay off.
set -u
cd "$(dirname "$0")"
export OMP_NUM_THREADS=12
while [ ! -f BTREF_DONE ]; do sleep 20; done      # wait for the restart pair to free the cores
rm -f BTREFS_DONE
echo "=== btref-scratch 0->100 from scratch, both arms parallel, 12 threads each $(date +%H:%M:%S)"
( env ATM_UBUD_BALANCE=1 ATM_CELL_ROT_DIAG=1                 ../cli/atm config_bs_ctl.xml > bs_ctl.log 2>&1
  echo "  bs_ctl exit $?  nan $(grep -ci 'nan' bs_ctl.log)  $(date +%H:%M:%S)" ) &
( env ATM_UBUD_BALANCE=1 ATM_CELL_ROT_DIAG=1 ATM_BUOY_TREF=1 ../cli/atm config_bs_on.xml  > bs_on.log  2>&1
  echo "  bs_on  exit $?  nan $(grep -ci 'nan' bs_on.log)   $(date +%H:%M:%S)" ) &
wait
echo "=== btref-scratch done $(date +%H:%M:%S)"
touch BTREFS_DONE
