#!/bin/bash
# BOTH-DIRECTIONS VERIFICATION OF THE FOUR DEFAULT FLIPS (2026-09-12).
#   ATM_CELLS_FROM_PSI            0   -> 1
#   ATM_TROPO_INDEX_FIX           0   -> 1
#   ATM_RADIAL_SHAPIRO_STRENGTH_VW (=STRENGTH) -> 0.0
#   ATM_V_MASSBAL_STRIDE          0   -> 1
# 1 THREAD, nm=20 FROM SCRATCH -- from scratch is mandatory: load_state() overwrites everything
# CELLS_FROM_PSI and TROPO_INDEX_FIX produce, so a restart cannot exercise either of them.
# 1 thread because this model is not bit-reproducible under OpenMP and a 24-thread pair could not
# tell a null flip from the threads.
#   A: new binary CLEAN     == old binary with the four SET      (the flip is what was measured)
#   B: new binary with four SET BACK == old binary CLEAN         (the revert restores the branch)
#   C: old+SET vs old+CLEAN -- the control; these MUST differ or A and B prove nothing.
#
# PARALLEL SINCE 2026-09-14. The four arms are independent 1-thread processes writing to four
# separate directories, so concurrency cannot change a byte and takes the wall clock from ~4T to
# ~T. See verify_lib.sh. DO NOT RE-SERIALISE.
set -u
cd "$(dirname "$0")"
. ./verify_lib.sh
NEW=../cli/atm; OLD=../cli/atm_preflip
ON="ATM_CELLS_FROM_PSI=1 ATM_TROPO_INDEX_FIX=1 ATM_RADIAL_SHAPIRO_STRENGTH_VW=0 ATM_V_MASSBAL_STRIDE=1"
OFF="ATM_CELLS_FROM_PSI=0 ATM_TROPO_INDEX_FIX=0 ATM_RADIAL_SHAPIRO_STRENGTH_VW=1.0 ATM_V_MASSBAL_STRIDE=0 ATM_VTK_STRIDE=1"
echo "=== launching 4 arms in parallel $(date +%H:%M:%S)"
arm va_old "$OLD" config_va_old.xml $ON
arm va_new "$NEW" config_va_new.xml
arm vb_old "$OLD" config_vb_old.xml
arm vb_new "$NEW" config_vb_new.xml $OFF
wait_arms
cmp_dirs    va_old va_new "A(flip)  new CLEAN == old+SET"
cmp_dirs    vb_old vb_new "B(revert) new+SETBACK == old CLEAN"
want_differ va_old vb_old "C(control) old+SET vs old CLEAN"
echo "=== verify done $(date +%H:%M:%S)"
