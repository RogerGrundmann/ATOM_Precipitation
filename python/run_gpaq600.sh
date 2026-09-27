#!/bin/bash
# Convection repairs scored on a WATER-CONSERVING branch (2026-09-26, POSTPONED -- launch by hand). 600 from scratch,
# cli/atm_rkp, 24 threads, gated on run_verify_vrks.sh (queue it after run_rksq100.sh, or run it alone).
#   gpaqb qh600 setup (ATM_WATER_CLOSURE=1, moisture filter off, T vert off) + the gpa_b convection stack
#         (ATM_MC_SGZ=1 ATM_MC_ENTR=1.0e-4 ATM_MC_QVD=2 ATM_MC_ALF1=5.44e-4 ATM_MC_BASE_SAT=2)
#         2026-09-27, at the user's instruction: BASE_SAT=2 ALONE, no GP_AREA -- gpa_b vs gpa_bc showed GP_AREA
#         makes convective rain 90x LARGER (8.8 -> 779) on top of (b). Renamed gpaq -> gpaqb for that reason.
# WHY: on the default branch the surface humidity is re-pinned, so precipitable water stays 30.4 mm and every mm of
# convective rain is refilled (gpa_bc: P - E 1187 mm/a, P_conv 779 on top of an unchanged ~1000 of other rain). With
# the closure, E is a flux (qh600: P/E 1.07) and convective rain must compete with stratiform rain for the same water.
# PRE-REGISTERED against qh600 (829 mm/a, r 0.455, sigma 2.23, bands 3044/115/91/2.9, land/ocean 779/849, P/E 1.07-1.11):
# (1) P/E stays ~1.1 -- the closure holds with convection active; (2) convective rain now REPLACES part of the other rain
# rather than adding to it: total precip stays near qh600's 829, not +779; (3) the open question: does that shift
# improve r / sigma / 0-15 / land? (4) exit 0, zero NaN through 155/357/483; MC_t truncation; max T location.
set -u; cd "$(dirname "$0")"; rm -f GPAQ600_DONE
until [ -f VRKS_VERIFY_DONE ]; do sleep 30; done
BAD=$(sed -n '1,/A  OFF BRANCH/p' run_verify_vrks.out | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "A  OFF BRANCH" run_verify_vrks.out; then
    echo "byte check did not pass -- NOT started"; cat run_verify_vrks.out; touch GPAQ600_DONE; exit 1; fi
mkdir output_gpaqb || { touch GPAQ600_DONE; exit 1; }
echo "gpaqb start $(date +%H:%M)"
env OMP_NUM_THREADS=${NT:-24} ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_MC_CAP_DIAG=1 \
    ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0 \
    ATM_MC_SGZ=1 ATM_MC_ENTR=1.0e-4 ATM_MC_QVD=2 ATM_MC_ALF1=5.44e-4 ATM_MC_BASE_SAT=2 \
    ${ATMBIN:-../cli/atm_rkp} config_gpaqb.xml > gpaqb.log 2>&1
echo "gpaqb exit $?  NaN $(grep -c 'NaN/Inf DETECTED' gpaqb.log)  $(date +%H:%M)"
touch GPAQ600_DONE
