#!/bin/bash
# OFF-BRANCH BYTE CHECK for HYD_VW_BOTTOM_ZG (new 2026-09-24, default 0), cli/hyd_zg (38fccc9a).
# 1 thread, nm = 20 from scratch, atm forcing seeded. Compared against the runs of run_verify_buoynd.sh:
#   bzg_new  hyd_zg unset                  -- MUST equal bnd_old (cli/hyd_bc2, HEAD)
#   bzg_bon  hyd_zg HYD_BUOY_CONSISTENT=1   -- MUST equal bnd_on  (cli/hyd_buoy, same knob)
#   bzg_zg   hyd_zg HYD_VW_BOTTOM_ZG=1      -- CONTROL: MUST differ from bzg_new
set -u; cd "$(dirname "$0")"; rm -f ZG_VERIFY_DONE
. ./verify_lib.sh
for t in new bon zg; do rm -rf output_bzg_$t; mkdir -p output_bzg_$t
  cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_bzg_$t/; done
arm bzg_new ../cli/hyd_zg config_bzg_new.xml
arm bzg_bon ../cli/hyd_zg config_bzg_bon.xml HYD_BUOY_CONSISTENT=1
arm bzg_zg  ../cli/hyd_zg config_bzg_zg.xml  HYD_VW_BOTTOM_ZG=1
wait_arms
cmp_dirs    bzg_new bnd_old "HYD A  OFF-BRANCH (hyd_zg unset vs HEAD)"
cmp_dirs    bzg_bon bnd_on  "HYD B  BUOY=1 unchanged by this build"
want_differ bzg_zg  bzg_new "HYD C  CONTROL ZG=1 vs unset -- MUST differ"
touch ZG_VERIFY_DONE
