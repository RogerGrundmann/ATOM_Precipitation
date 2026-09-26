#!/bin/bash
# The storm-track warming, next suspect (2026-09-26, POSTPONED -- launch by hand). 100 from scratch, cli/atm_rkp,
# 24 threads, gated on run_verify_vrks.sh.
#   rksq_on  qh_on setup (closure + moisture filter off + T vert off) + ATM_RK_SCALAR_SYNC=1 (moisture only)
# Against qh_on (sync 2) and qh_off (closure off). qh_on's 35-65 column is +0.65 K warmer than qh_off at every height,
# with the same vapour, half the cloud and half the rain; qth_on exonerated the T filter's horizontal passes.
# PRE-REGISTERED: if the pre-RK4 SaturationAdjustment latent heat (kept by mode 2) is the warming, the 35-65 column
# falls back toward qh_off (-0.65 K), cloud water toward 28.8 g/m2, 15-35 / 35-65 precip toward 217 / 192.
# CONFOUND, STATED BEFORE THE RUN: mode 1 keeps the condensate but discards its latent heat -- "physically
# inconsistent, a probe" -- so SOME cooling is expected by construction. A cooling does not prove the mechanism by
# itself; NO cooling would refute it. Read with t_z.py / band_moist.py on output_rksq_on/atm_restart_0Ma_100.bin.
set -u; cd "$(dirname "$0")"; rm -f RKSQ100_DONE
until [ -f VRKS_VERIFY_DONE ]; do sleep 30; done
BAD=$(sed -n '1,/B  CLOSURE/p' run_verify_vrks.out | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "CONTROL.*PASS" run_verify_vrks.out; then
    echo "byte check did not pass -- NOT started"; cat run_verify_vrks.out; touch RKSQ100_DONE; exit 1; fi
mkdir output_rksq_on || { touch RKSQ100_DONE; exit 1; }
echo "rksq_on start $(date +%H:%M)"
env OMP_NUM_THREADS=24 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_MC_CAP_DIAG=1 \
    ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0 ATM_RK_SCALAR_SYNC=1 \
    ../cli/atm_rkp config_rksq_on.xml > rksq_on.log 2>&1
echo "rksq_on exit $?  NaN $(grep -c 'NaN/Inf DETECTED' rksq_on.log)  $(date +%H:%M)"
touch RKSQ100_DONE
