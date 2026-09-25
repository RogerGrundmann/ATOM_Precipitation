#!/bin/bash
# OFF-BRANCH BYTE CHECK for ATM_DAMP_Q_VERT (default 1 = shipped), 2026-09-25. 1 thread, nm = 20 from scratch.
# old = cli/atm_wcr (closure default off), new = cli/atm_dqv (+ the knob).
set -u; cd "$(dirname "$0")"; rm -f DQV_VERIFY_DONE
. ./verify_lib.sh
for t in old new ctl; do mkdir output_dqv_$t || exit 1; done
arm dqv_old ../cli/atm_wcr config_dqv_old.xml
arm dqv_new ../cli/atm_dqv config_dqv_new.xml
arm dqv_ctl ../cli/atm_dqv config_dqv_ctl.xml ATM_DAMP_Q_VERT=0
wait_arms
cmp_dirs    dqv_new dqv_old "A  OFF BRANCH (new clean vs old)"
want_differ dqv_ctl dqv_new "C  CONTROL -- MUST differ"
diff <(sed "s#output_dqv_new#X#g" output_dqv_new/RUN_CONFIG.txt) <(sed "s#output_dqv_old#X#g" output_dqv_old/RUN_CONFIG.txt) | grep -o "DAMP_Q_VERT=[^ ]*" | tr '\n' ' '; echo
touch DQV_VERIFY_DONE
