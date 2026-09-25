#!/bin/bash
# The closure's POLAR engine, located (2026-09-25): qvt_on repeated with ATM_CWB_BANDS=1. 100 from scratch,
# cli/atm_cwbb, 6 threads. ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_T_VERT=0 -- the arm where 65-90 still
# climbed 28 -> 669 mm/a. Must reproduce qvt_on (the instrument is print-only). Read the 65-90 column of
# [CWB-BANDS]: which stage puts water into the polar band, against the closure-off shape (bis_wc).
set -u; cd "$(dirname "$0")"; rm -f QVTB100_DONE
mkdir output_qvtb_on || { touch QVTB100_DONE; exit 1; }
echo "start $(date +%H:%M)"
env OMP_NUM_THREADS=6 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_MC_CAP_DIAG=1 ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_T_VERT=0 \
    ../cli/atm_cwbb config_qvtb_on.xml > qvtb_on.log 2>&1
echo "qvtb_on exit $?  NaN $(grep -c 'NaN/Inf DETECTED' qvtb_on.log)  $(date +%H:%M)"
touch QVTB100_DONE
