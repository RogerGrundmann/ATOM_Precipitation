#!/bin/bash
# Arms for items 1, 2, 6 (2026-09-26, POSTPONED -- launch by hand). 100 from scratch, cli/atm_fx4, 6 arms x 4 threads
# concurrent, gated on run_verify_vfx4.sh. ATM_CWB_DIAG + ATM_CWB_BANDS on every arm.
#   fx4_ctl     default                                  | fx4_tsf  default + ATM_TURB_SIN_FLOOR=1
#   fx4_sqcoff  qh_on setup (closure on, moisture filter off)   | fx4_sqcon  same + ATM_SEAM_Q_CONSERVE=1
#   fx4_oc1     OneCat (CategoryIceScheme 1)             | fx4_oc2  OneCat + ATM_ONECAT_CLOUD_LIMIT=1
# PRE-REGISTERED: (2) tsf: only rows poleward of ~75 deg move (turbulence production/nue there); precip bands
#   0-65 identical to the digit, 65-90 may move; exit 0, zero NaN. (6) sqcon: CWB "BC:phi(seam)" goes to ~0
#   (the self-check; clip residual only) where sqcoff shows the leak (+8.2e3 into 15-35 in qh_on); climate
#   otherwise ~null. (1) oc2: CWB RungeKutta clip bucket for cloud shrinks against oc1; OneCat precip moves.
set -u; cd "$(dirname "$0")"; rm -f FX4_100_DONE
until [ -f VFX4_VERIFY_DONE ]; do sleep 30; done
BAD=$(sed -n '1,/A  OFF BRANCH/p' run_verify_vfx4.out | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "A  OFF BRANCH" run_verify_vfx4.out || [ "$(grep -c 'CONTROL.*PASS' run_verify_vfx4.out)" != 3 ]; then
    echo "byte check did not pass -- NOT started"; cat run_verify_vfx4.out; touch FX4_100_DONE; exit 1; fi
for t in ctl tsf sqcoff sqcon oc1 oc2; do mkdir output_fx4_$t || { touch FX4_100_DONE; exit 1; }; done
QH="ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0"
echo "start $(date +%H:%M)"
run(){ ( env OMP_NUM_THREADS=4 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_MC_CAP_DIAG=1 $2 ../cli/atm_fx4 config_$1.xml > $1.log 2>&1
         echo "$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $1.log)  $(date +%H:%M)" ) & }
run fx4_ctl    ""
run fx4_tsf    "ATM_TURB_SIN_FLOOR=1"
run fx4_sqcoff "$QH"
run fx4_sqcon  "$QH ATM_SEAM_Q_CONSERVE=1"
run fx4_oc1    ""
run fx4_oc2    "ATM_ONECAT_CLOUD_LIMIT=1"
wait
touch FX4_100_DONE
