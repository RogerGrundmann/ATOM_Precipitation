#!/bin/bash
# OFF-BRANCH BYTE CHECK for ATM_DAMP_Q_HORIZ (default 1 = shipped), 2026-09-26. 1 thread, nm = 20 from scratch.
# old = cli/atm_predqh (b46f036), new = cli/atm_dqh (+ the knob). Control: the moisture filter fully off.
set -u; cd "$(dirname "$0")"; rm -f DQH_VERIFY_DONE
. ./verify_lib.sh
for t in old new ctl; do mkdir output_dqh_$t || exit 1; done
arm dqh_old ../cli/atm_predqh config_dqh_old.xml
arm dqh_new ../cli/atm_dqh    config_dqh_new.xml
arm dqh_ctl ../cli/atm_dqh    config_dqh_ctl.xml ATM_DAMP_Q_HORIZ=0
wait_arms
cmp_dirs    dqh_new dqh_old "A  OFF BRANCH (new clean vs old)"
want_differ dqh_ctl dqh_new "C  CONTROL -- MUST differ"
diff <(sed "s#output_dqh_new#X#g" output_dqh_new/RUN_CONFIG.txt) <(sed "s#output_dqh_old#X#g" output_dqh_old/RUN_CONFIG.txt) | grep -o "DAMP_Q_HORIZ=[^ ]*" | tr '\n' ' '; echo
touch DQH_VERIFY_DONE
