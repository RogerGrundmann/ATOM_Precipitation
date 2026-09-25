#!/bin/bash
# OFF-BRANCH BYTE CHECK for ATM_DAMP_T_VERT (default 1 = shipped), 2026-09-25. 1 thread, nm = 20 from scratch.
# old = cli/atm_dqv (1fe8475), new = cli/atm_dtv (+ the knob).
set -u; cd "$(dirname "$0")"; rm -f DTV_VERIFY_DONE
. ./verify_lib.sh
for t in old new ctl; do mkdir output_dtv_$t || exit 1; done
arm dtv_old ../cli/atm_dqv config_dtv_old.xml
arm dtv_new ../cli/atm_dtv config_dtv_new.xml
arm dtv_ctl ../cli/atm_dtv config_dtv_ctl.xml ATM_DAMP_T_VERT=0
wait_arms
cmp_dirs    dtv_new dtv_old "A  OFF BRANCH (new clean vs old)"
want_differ dtv_ctl dtv_new "C  CONTROL -- MUST differ"
diff <(sed "s#output_dtv_new#X#g" output_dtv_new/RUN_CONFIG.txt) <(sed "s#output_dtv_old#X#g" output_dtv_old/RUN_CONFIG.txt) | grep -o "DAMP_T_VERT=[^ ]*" | tr '\n' ' '; echo
touch DTV_VERIFY_DONE
