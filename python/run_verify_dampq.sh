#!/bin/bash
# ATM_DAMP_Q_MASS (B+4, new 2026-09-24, default 0). SHORT CHECKS ONLY -- long arms postponed by the user.
# (1) off-branch byte check, 1 thread, nm = 20 from scratch:
#   bdq_old  cli/atm_sd (ee2ddaa4-equivalent HEAD 14b7cad atmosphere), clean
#   bdq_new  cli/atm_dq, unset -- MUST equal bdq_old except RUN_CONFIG.txt (banner token)
#   bdq_on   cli/atm_dq, =1    -- CONTROL: MUST differ
# (2) conservation: 600 -> 620 restarts from output_hsf_ctl's checkpoint with ATM_CWB_DIAG=1,
#     knob 0 vs 1, 3 threads each. The damp_wiggles(q) bucket must fall from ~2.5e+06 mm/a to ~0.
set -u; cd "$(dirname "$0")"; rm -f DQ_VERIFY_DONE
. ./verify_lib.sh
for t in old new on; do rm -rf output_bdq_$t; mkdir -p output_bdq_$t; done
for t in c0 c1; do mkdir -p output_dqc_$t; cp -n output_hsf_ctl/atm_restart_0Ma_600.bin output_dqc_$t/
  sed -e "s#output_hsf1200_ctl/#output_dqc_$t/#" -e 's#<nm>1200</nm>#<nm>620</nm>#' \
      -e 's#<checkpoint_save_iter>1200<#<checkpoint_save_iter>-1<#' config_hsf1200_ctl.xml > config_dqc_$t.xml; done
arm bdq_old ../cli/atm_sd config_bdq_old.xml
arm bdq_new ../cli/atm_dq config_bdq_new.xml
arm bdq_on  ../cli/atm_dq config_bdq_on.xml ATM_DAMP_Q_MASS=1
( env OMP_NUM_THREADS=3 ATM_CWB_DIAG=1                   ../cli/atm_dq config_dqc_c0.xml > dqc_c0.log 2>&1 ) &
( env OMP_NUM_THREADS=3 ATM_CWB_DIAG=1 ATM_DAMP_Q_MASS=1 ../cli/atm_dq config_dqc_c1.xml > dqc_c1.log 2>&1 ) &
wait_arms; wait
cmp_dirs    bdq_new bdq_old "ATM A  OFF-BRANCH (atm_dq unset vs HEAD)"
diff <(sed 's#output_bdq_new#X#g' output_bdq_new/RUN_CONFIG.txt) <(sed 's#output_bdq_old#X#g' output_bdq_old/RUN_CONFIG.txt) | grep -o "DAMP_Q_MASS=[^ ]*"
want_differ bdq_on  bdq_new "ATM C  CONTROL =1 vs unset -- MUST differ"
for t in c0 c1; do echo "== dqc_$t  exit NaN=$(grep -c 'NaN/Inf DETECTED' dqc_$t.log)"; grep -E "damp_wiggles\(q\)|unattributed|NET|Precip mean" dqc_$t.log | tail -5; done
touch DQ_VERIFY_DONE
