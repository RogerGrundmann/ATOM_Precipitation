#!/bin/bash
# OFF-BRANCH BYTE CHECK for the RK_SCALAR_SYNC precedence fix (2026-09-26). 1 thread, nm = 20 from scratch.
# old = cli/atm_fx4 (476196e source), new = cli/atm_rkp. LAUNCH BY HAND (postponed).
#   A  clean: new == old
#   B  ATM_WATER_CLOSURE=1 alone: new == old (the closure still forces mode 2 when RK_SCALAR_SYNC is unset)
#   C  closure + ATM_RK_SCALAR_SYNC=1 on new MUST differ from closure alone (the override now takes effect)
set -u; cd "$(dirname "$0")"; rm -f VRKS_VERIFY_DONE
. ./verify_lib.sh
for t in old new wco wcn ctl; do mkdir output_vrks_$t || { touch VRKS_VERIFY_DONE; exit 1; }; done
arm vrks_old ../cli/atm_fx4 config_vrks_old.xml
arm vrks_new ../cli/atm_rkp config_vrks_new.xml
arm vrks_wco ../cli/atm_fx4 config_vrks_wco.xml ATM_WATER_CLOSURE=1
arm vrks_wcn ../cli/atm_rkp config_vrks_wcn.xml ATM_WATER_CLOSURE=1
arm vrks_ctl ../cli/atm_rkp config_vrks_ctl.xml ATM_WATER_CLOSURE=1 ATM_RK_SCALAR_SYNC=1
wait_arms
cmp_dirs    vrks_new vrks_old "A  OFF BRANCH (new clean vs old)"
cmp_dirs    vrks_wcn vrks_wco "B  CLOSURE ALONE unchanged (new vs old)"
want_differ vrks_ctl vrks_wcn "C  CONTROL closure + RK_SCALAR_SYNC=1 -- MUST differ"
touch VRKS_VERIFY_DONE
