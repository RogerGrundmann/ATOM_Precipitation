#!/bin/bash
# BOTH-DIRECTIONS VERIFICATION OF THE ATM_BUOY_CONSISTENT REVERT (2026-09-14).
#   ATM_BUOY_CONSISTENT  1 (default 2026-09-09 .. 2026-09-13)  ->  0 (default again)
# 1 THREAD, nm=20 FROM SCRATCH. 1 thread because this model is not bit-reproducible under OpenMP
# and a 24-thread pair could not tell a null flip from the threads. From scratch because a restart
# would start both arms from a field one of them could not have produced.
#   A: new binary CLEAN        == old binary with ATM_BUOY_CONSISTENT=0   (the revert)
#   B: new binary with =1 SET  == old binary CLEAN                        (the flipped branch)
#   C: old+=0 vs old CLEAN -- the control; these MUST differ or A and B prove nothing.
#      buoyancy_ramp = min(1, iter/300), so at nm=20 the term runs at 0.3-6.7 % strength -- small
#      but non-zero from iteration 1, against a branch difference of 5.0e5. It fires.
# PARALLEL BY DEFAULT -- see verify_lib.sh. DO NOT RE-SERIALISE.
set -u
cd "$(dirname "$0")"
. ./verify_lib.sh
NEW=../cli/atm; OLD=../cli/atm_bcflip_pre
echo "=== launching 4 arms in parallel $(date +%H:%M:%S)"
arm ba_old "$OLD" config_ba_old.xml ATM_BUOY_CONSISTENT=0
arm ba_new "$NEW" config_ba_new.xml
arm bb_old "$OLD" config_bb_old.xml
arm bb_new "$NEW" config_bb_new.xml ATM_BUOY_CONSISTENT=1
wait_arms
cmp_dirs    ba_old ba_new "A(revert)  new CLEAN == old+=0"
cmp_dirs    bb_old bb_new "B(flipped) new+=1 == old CLEAN"
want_differ ba_old bb_old "C(control) old+=0 vs old CLEAN"
echo "=== verify done $(date +%H:%M:%S)"
touch VERIFY_BUOY_DONE
