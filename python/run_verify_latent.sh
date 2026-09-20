#!/bin/bash
# ATM_SATADJ_FREEZE_LATENT both-directions byte check. Waits for the FADE check to clear first,
# because rebuilding cli/atm while that check's arms are executing it is a build hazard.
#   vl_new  new binary, unset              -- must equal vl_old
#   vl_old  cli/atm_sf, the pre-knob binary
#   vl_on   new binary, =1                 -- want_differ CONTROL. If this reads VACUOUS the knob
#                                             found no cell below t_00 holding liquid in 20
#                                             iterations, and the check has to be run longer.
set -u; cd "$(dirname "$0")"
while [ ! -f VFADE_DONE ]; do sleep 30; done
echo "=== fade check cleared, rebuilding $(date +%H:%M:%S)"
( cd .. && make atm ) > build_latent.log 2>&1 || { echo "BUILD FAILED"; tail -5 build_latent.log; exit 1; }
echo "=== build ok $(date +%H:%M:%S)"
. ./verify_lib.sh
rm -rf output_vl_new output_vl_old output_vl_on; mkdir -p output_vl_new output_vl_old output_vl_on
arm vl_new ../cli/atm    config_vl_new.xml
arm vl_old ../cli/atm_sf config_vl_old.xml
arm vl_on  ../cli/atm    config_vl_on.xml  ATM_SATADJ_FREEZE_LATENT=1
wait_arms
cmp_dirs    vl_new vl_old "A  off-branch (new unset vs pre-knob cli/atm_sf)"
want_differ vl_on  vl_new "C  CONTROL freeze-latent=1 vs off -- MUST differ"
touch VLATENT_DONE
