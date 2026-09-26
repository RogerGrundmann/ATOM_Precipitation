#!/bin/bash
# OFF-BRANCH BYTE CHECK for HYD_METRIC_SIN_FLOOR (item 3, default 0.4 = shipped), 2026-09-26. 1 thread, nm = 20
# from scratch, atm forcing seeded. old = cli/hyd_hsp (b6f457e), new = cli/hyd_fx4. LAUNCH BY HAND (postponed).
set -u; cd "$(dirname "$0")"; rm -f VOSF_VERIFY_DONE
. ./verify_lib.sh
for t in old new ctl; do mkdir output_vosf_$t || { touch VOSF_VERIFY_DONE; exit 1; }
  cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_vosf_$t/; done
arm vosf_old ../cli/hyd_hsp config_vosf_old.xml
arm vosf_new ../cli/hyd_fx4 config_vosf_new.xml
arm vosf_ctl ../cli/hyd_fx4 config_vosf_ctl.xml HYD_METRIC_SIN_FLOOR=0.26
wait_arms
cmp_dirs    vosf_new vosf_old "HYD A  OFF-BRANCH (new unset vs b6f457e)"
want_differ vosf_ctl vosf_new "HYD C  CONTROL METRIC_SIN_FLOOR=0.26 -- MUST differ"
touch VOSF_VERIFY_DONE
