#!/bin/bash
# STEP A, long arm (2026-09-26): closure ON with the moisture filter fully OFF, 600 from scratch, cli/atm_dqh,
# 24 threads, queued behind run_alf600.sh (ALF600_DONE). Reference: sgzb_ctl (corrected default, closure off,
# 600 from scratch: 987 mm/a, drift -0.01/iter, r 0.471, sigma 2.37, bands 3388/211/205/24, land/ocean 854/1039).
#   qh600  ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0   (= qh_on to 600)
# PRE-REGISTERED (qh_on at 100: 818 mm/a rising 3.26/iter, 15-35 116, 35-65 92, 65-90 2.9, land 774, P/E 787/757):
# (1) the question: do 15-35 / 35-65 recover toward sgzb_ctl's 211 / 205, or stay starved at ~half? No meridional
#     moisture supply is left but the resolved transport, so STARVED is the expectation.
# (2) the global climb must LEVEL OFF (drift 400-600 near 0); a sustained climb means a third engine.
# (3) 65-90 stays below ~30 (no polar flooding). (4) CWB damp_wiggles(q) exactly 0; P - E stays small.
# (5) exit 0, zero NaN through 155/357/483; no grid-scale moisture mode (max cloud water location/value).
set -u; cd "$(dirname "$0")"; rm -f QH600_DONE
until [ -f ALF600_DONE ]; do sleep 60; done
mkdir output_qh600 || { touch QH600_DONE; exit 1; }
echo "qh600 start $(date +%H:%M)"
env OMP_NUM_THREADS=24 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_MC_CAP_DIAG=1 \
    ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0 \
    ../cli/atm_dqh config_qh600.xml > qh600.log 2>&1
echo "qh600 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' qh600.log)  $(date +%H:%M)"
touch QH600_DONE
