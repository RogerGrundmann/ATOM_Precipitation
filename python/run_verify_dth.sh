#!/bin/bash
# OFF-BRANCH BYTE CHECK for ATM_DAMP_T_HORIZ (default 1 = shipped), 2026-09-26. 1 thread, nm = 20 from scratch.
# old = cli/atm_alf (425680f), new = cli/atm_dth (+ the knob). Queued behind run_qh600.sh (QH600_DONE).
set -u; cd "$(dirname "$0")"; rm -f DTH_VERIFY_DONE
until [ -f QH600_DONE ]; do sleep 60; done
. ./verify_lib.sh
for t in old new ctl; do mkdir output_dth_$t || { touch DTH_VERIFY_DONE; exit 1; }; done
arm dth_old ../cli/atm_alf config_dth_old.xml
arm dth_new ../cli/atm_dth config_dth_new.xml
arm dth_ctl ../cli/atm_dth config_dth_ctl.xml ATM_DAMP_T_HORIZ=0
wait_arms
cmp_dirs    dth_new dth_old "A  OFF BRANCH (new clean vs old)"
want_differ dth_ctl dth_new "C  CONTROL -- MUST differ"
diff <(sed "s#output_dth_new#X#g" output_dth_new/RUN_CONFIG.txt) <(sed "s#output_dth_old#X#g" output_dth_old/RUN_CONFIG.txt) | grep -o "DAMP_T_HORIZ=[^ ]*" | tr '\n' ' '; echo
touch DTH_VERIFY_DONE
