#!/bin/bash
# STEP A (2026-09-26): the moisture Shapiro filter FULLY OFF under the closure. 100 from scratch, cli/atm_dqh,
# 8 threads each, both arms concurrent, gated on run_verify_dqh.sh.
#   qh_on   ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0
#           -- one variable against qvtb_on (closure on, both vertical passes off, horizontal ON)
#   qh_off  ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0   -- closure off: the filter's effect on the shipped branch (vs bis_wc)
# PRE-REGISTERED: (1) qh_on's damp_wiggles(q) CWB bucket is EXACTLY 0 (self-check). (2) the 65-90 band stays
# flat near the closure-off 12-16 mm/a instead of qvtb_on's 29 -> 675. (3) the global climb over 40-100 falls to
# the closure-off ~3.6 mm/a/iter or below (qv_on 5.44). (4) exit 0, no NaN: the filter was installed against 2dt
# moisture modes -- watch max cloud water / max vapour and their location for a grid-scale runaway.
# (5) qh_off is a near-null against bis_wc (RK4 discards the filter on that branch).
set -u; cd "$(dirname "$0")"; rm -f QH100_DONE
until [ -f DQH_VERIFY_DONE ]; do sleep 20; done
BAD=$(sed -n '1,/A  OFF BRANCH/p' run_verify_dqh.out | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "A  OFF BRANCH" run_verify_dqh.out || ! grep -q "CONTROL.*PASS" run_verify_dqh.out; then
    echo "byte check did not pass -- NOT started"; cat run_verify_dqh.out; touch QH100_DONE; exit 1; fi
for t in on off; do mkdir output_qh_$t || { touch QH100_DONE; exit 1; }; done
echo "start $(date +%H:%M)"
( env OMP_NUM_THREADS=8 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_MC_CAP_DIAG=1 ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0 \
    ../cli/atm_dqh config_qh_on.xml > qh_on.log 2>&1; echo "qh_on exit $?  NaN $(grep -c 'NaN/Inf DETECTED' qh_on.log)  $(date +%H:%M)" ) &
( env OMP_NUM_THREADS=8 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_MC_CAP_DIAG=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 \
    ../cli/atm_dqh config_qh_off.xml > qh_off.log 2>&1; echo "qh_off exit $?  NaN $(grep -c 'NaN/Inf DETECTED' qh_off.log)  $(date +%H:%M)" ) &
wait
touch QH100_DONE
