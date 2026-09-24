#!/bin/bash
# ATM_WATER_CLOSURE (B+3 flux + RK4 sync mode 2, new 2026-09-24, default 0). SHORT CHECKS ONLY --
# longer runs postponed by the user.
# (1) byte check, 1 thread, nm = 20 from scratch: bwc_new (atm_wc unset) MUST equal bdq_new (atm_dq
#     unset; that one is checked against HEAD by run_verify_dampq.sh) except RUN_CONFIG.txt;
#     bwc_on (=1) MUST differ.
# (2) 600 -> 620 restart from output_hsf_ctl with ATM_CWB_DIAG=1, knob on (dwc_c1); control is
#     dqc_c0 (run_verify_dampq.sh, off branch). Self-checks: leapfrog_reset = 0 exactly; the
#     evaporation bucket ~ the printed E (hundreds of mm/a), not ~5e+06; zero NaN.
set -u; cd "$(dirname "$0")"; rm -f WC_VERIFY_DONE
until [ -f DQ_VERIFY_DONE ]; do sleep 30; done
. ./verify_lib.sh
for t in new on; do rm -rf output_bwc_$t; mkdir -p output_bwc_$t; done
arm bwc_new ../cli/atm_wc config_bwc_new.xml
arm bwc_on  ../cli/atm_wc config_bwc_on.xml ATM_WATER_CLOSURE=1
( env OMP_NUM_THREADS=3 ATM_CWB_DIAG=1 ATM_WATER_CLOSURE=1 ../cli/atm_wc config_dwc_c1.xml > dwc_c1.log 2>&1 ) &
wait_arms; wait
cmp_dirs    bwc_new bdq_new "ATM A  OFF-BRANCH (atm_wc unset vs atm_dq unset)"
want_differ bwc_on  bwc_new "ATM C  CONTROL =1 vs unset -- MUST differ"
for t in dqc_c0 dwc_c1; do echo "== $t  NaN=$(grep -c 'NaN/Inf DETECTED' $t.log)"
  grep -E "evaporation|RungeKutta|leapfrog_reset|damp_wiggles\(q\)|unattributed|NET|water budget closure|Precip mean|WATER CLOSURE" $t.log | tail -12; done
touch WC_VERIFY_DONE
