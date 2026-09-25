#!/bin/bash
# BOTH vertical filter passes off under the closure (2026-09-25). 100 from scratch, cli/atm_dtv, 6 threads,
# queued behind run_qv100.sh and gated on run_verify_dtv.sh. Compare bis_ctl (closure on, both passes on),
# qv_on (closure on, moisture pass off) and bis_wc (closure off).
#   qvt_on  ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_T_VERT=0
# PRE-REGISTERED: climb over 40-100 at or below qv_on's; if qv_on already sits at the closure-off level
# (~3.6 mm/a/iter) this arm is expected to be a null against it. Watch the Pamir surface t (the reason
# the t filter was added: a 2dt mode that ran to 53 K): max temperature / Pamir t in the log, exit 0, no NaN.
set -u; cd "$(dirname "$0")"; rm -f QVT100_DONE
until [ -f QV100_DONE ] && [ -f DTV_VERIFY_DONE ]; do sleep 30; done
BAD=$(sed -n '1,/A  OFF BRANCH/p' run_verify_dtv.out | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "A  OFF BRANCH" run_verify_dtv.out || ! grep -q "CONTROL.*PASS" run_verify_dtv.out; then
    echo "byte check did not pass -- NOT started"; cat run_verify_dtv.out; touch QVT100_DONE; exit 1; fi
mkdir output_qvt_on || { touch QVT100_DONE; exit 1; }
echo "start $(date +%H:%M)"
env OMP_NUM_THREADS=6 ATM_CWB_DIAG=1 ATM_MC_CAP_DIAG=1 ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_T_VERT=0 \
    ../cli/atm_dtv config_qvt_on.xml > qvt_on.log 2>&1
echo "qvt_on exit $?  NaN $(grep -c 'NaN/Inf DETECTED' qvt_on.log)  $(date +%H:%M)"
touch QVT100_DONE
