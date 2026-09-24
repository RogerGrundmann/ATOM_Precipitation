#!/bin/bash
# OFF-BRANCH BYTE CHECK for HYD_BUOY_CONSISTENT (new 2026-09-24, default 0.0). 1 thread, nm = 20 from
# scratch, atm forcing output_twctl/..._Transfer_Atm_600.vwtp seeded into every arm.
#   bnd_old  cli/hyd_bc2 (63d3c910, HEAD hydrosphere, knob absent)
#   bnd_new  cli/hyd_buoy (knob unset)          -- MUST equal bnd_old (RUN_CONFIG/banner aside)
#   bnd_on   cli/hyd_buoy, HYD_BUOY_CONSISTENT=1 -- CONTROL: MUST differ (the knob is connected)
set -u; cd "$(dirname "$0")"; rm -f BUOYND_VERIFY_DONE
. ./verify_lib.sh
for t in old new on; do rm -rf output_bnd_$t; mkdir -p output_bnd_$t
  cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_bnd_$t/; done
arm bnd_old ../cli/hyd_bc2  config_bnd_old.xml
arm bnd_new ../cli/hyd_buoy config_bnd_new.xml
arm bnd_on  ../cli/hyd_buoy config_bnd_on.xml HYD_BUOY_CONSISTENT=1
wait_arms
cmp_dirs    bnd_new bnd_old "HYD A  OFF-BRANCH (new unset vs HEAD)"
want_differ bnd_on  bnd_new "HYD C  CONTROL =1 vs unset -- MUST differ"
grep -h "buoyancy coefficient" bnd_new.log | head -1
grep -c "NaN/Inf DETECTED" bnd_on.log
touch BUOYND_VERIFY_DONE
