#!/bin/bash
# BOTH-DIRECTIONS BYTE CHECK for the 2026-09-23 flip of ATM_RAD_TOPO 0 -> 1 (default).
#   vrt_old  cli/atm_prert (pre-flip, rebuilt 2026-09-24 at 81b1251, md5 e02f5d29), clean  = shipped branch
#   vrt_off  cli/atm_rt (post-flip, rebuilt 2026-09-24, md5 4305c587), ATM_RAD_TOPO=0 -- MUST equal vrt_old
#   vrt_new  cli/atm_rt, clean                              = the new default
#   vrt_on   cli/atm_prert, ATM_RAD_TOPO=1                  -- MUST equal vrt_new
#   CONTROL  vrt_new vs vrt_off MUST DIFFER (MLR runs at setup and feeds t_eq via the CO2
#            perturbation, so it should differ from iteration 0).
# RUN_CONFIG.txt differs by banner token and path only. 1 thread, nm = 20 from scratch.
set -u; cd "$(dirname "$0")"; rm -f RADTOPO_VERIFY_DONE
. ./verify_lib.sh
rm -rf output_vrt_old output_vrt_off output_vrt_new output_vrt_on
mkdir -p output_vrt_old output_vrt_off output_vrt_new output_vrt_on
arm vrt_old ../cli/atm_prert config_vrt_old.xml
arm vrt_off ../cli/atm_rt    config_vrt_off.xml ATM_RAD_TOPO=0
arm vrt_new ../cli/atm_rt    config_vrt_new.xml
arm vrt_on  ../cli/atm_prert config_vrt_on.xml  ATM_RAD_TOPO=1
wait_arms
cmp_dirs    vrt_off vrt_old "A  RESTORE direction (new binary =0  vs  pre-flip clean)"
cmp_dirs    vrt_on  vrt_new "B  FLIP direction    (pre-flip =1    vs  new binary clean)"
want_differ vrt_new vrt_off "C  CONTROL new default vs off -- MUST differ"
touch RADTOPO_VERIFY_DONE
