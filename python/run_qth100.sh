#!/bin/bash
# The storm-track temperature test (2026-09-26): qh_on + the TEMPERATURE filter's horizontal passes off.
# 100 from scratch, cli/atm_dth, 24 threads, queued behind run_verify_dth.sh (itself behind qh600).
#   qth_on  ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0 ATM_DAMP_T_HORIZ=0
#           -- one variable against qh_on (8 thr; thread noise is ~1e-4, far below the effect sought)
# PRE-REGISTERED (qh_on at 100: 35-65 column +0.65 K vs qh_off at every height, LWP 15.8 vs 28.8 g/m2,
# precip 15-35/35-65 116/92 vs 217/192): if the horizontal t pass is the heat pump, the 35-65 column cools
# back toward qh_off (band_moist/t_z scratch scripts on output_qth_on/atm_restart_0Ma_100.bin), LWP and the two
# extratropical bands recover toward the closure-off values. If not, the next suspect is pre-RK4 SatAdj latent
# heat kept by the sync. STABILITY: this filter was added against a Pamir 2dt surface-t mode that ran to 53 K --
# watch max temperature and its location (Pamir ~38N 73E), exit 0, zero NaN.
set -u; cd "$(dirname "$0")"; rm -f QTH100_DONE
until [ -f DTH_VERIFY_DONE ]; do sleep 30; done
BAD=$(sed -n '1,/A  OFF BRANCH/p' run_verify_dth.out | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "A  OFF BRANCH" run_verify_dth.out || ! grep -q "CONTROL.*PASS" run_verify_dth.out; then
    echo "byte check did not pass -- NOT started"; cat run_verify_dth.out; touch QTH100_DONE; exit 1; fi
mkdir output_qth_on || { touch QTH100_DONE; exit 1; }
echo "qth_on start $(date +%H:%M)"
env OMP_NUM_THREADS=24 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_MC_CAP_DIAG=1 \
    ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0 ATM_DAMP_T_HORIZ=0 \
    ../cli/atm_dth config_qth_on.xml > qth_on.log 2>&1
echo "qth_on exit $?  NaN $(grep -c 'NaN/Inf DETECTED' qth_on.log)  $(date +%H:%M)"
touch QTH100_DONE
