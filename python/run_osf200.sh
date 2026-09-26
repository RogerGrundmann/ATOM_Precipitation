#!/bin/bash
# Item 3: the ocean polar metric floor, swept (2026-09-26, POSTPONED -- launch by hand). 1000 -> 1200 restarts from
# output_oc_ctl (shipped metric, correct salinity), cli/hyd_fx4, 3 arms x 8 threads, gated on run_verify_vosf.sh.
#   osf_40  0.4 (shipped, ~66 deg)   osf_26  0.26 (~75 deg, the atmosphere's value)   osf_10  0.10 (~84 deg)
# PRE-REGISTERED: the atmosphere's sweep was a null outside the rows where the floor binds (the polar filter did
# the work there). Expect the same: equatorward of 66 deg identical; poleward, max|v|/|w| and KE move by a few %.
# A blow-up at 0.10 (explicit 1/sin^2 diffusion) is the failure to look for: exit status, NaN, max|u|.
set -u; cd "$(dirname "$0")"; rm -f OSF200_DONE
until [ -f VOSF_VERIFY_DONE ]; do sleep 60; done
if ! grep -q "OFF-BRANCH.*PASS" run_verify_vosf.out || ! grep -q "CONTROL.*PASS" run_verify_vosf.out; then
    echo "byte check did not pass -- NOT started"; cat run_verify_vosf.out; touch OSF200_DONE; exit 1; fi
for t in 40 26 10; do mkdir output_osf_$t || { touch OSF200_DONE; exit 1; }
  cp output_oc_ctl/hyd_restart_0Ma_1000.bin output_osf_$t/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_osf_$t/; done
echo "start $(date +%H:%M)"
run(){ ( env OMP_NUM_THREADS=8 HYD_METRIC_SIN_FLOOR=$2 ../cli/hyd_fx4 config_osf_$1.xml > osf_$1.log 2>&1
         echo "osf_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' osf_$1.log)  $(date +%H:%M)" ) & }
run 40 0.4; run 26 0.26; run 10 0.10
wait
touch OSF200_DONE
