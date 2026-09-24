#!/bin/bash
# OFF-BRANCH BYTE CHECK for HYD_DEEP_DRAG (B.7, new 2026-09-24, default 0). 1 thread, nm = 20 from scratch,
# atm forcing seeded. Longer runs postponed at the user's instruction.
#   bdd_new  cli/hyd_dd unset        -- MUST equal bzg_new (cli/hyd_zg unset, HEAD ocean)
#   bdd_on   cli/hyd_dd, =1 (1 day)  -- CONTROL: MUST differ
set -u; cd "$(dirname "$0")"; rm -f DD_VERIFY_DONE
. ./verify_lib.sh
for t in new on; do rm -rf output_bdd_$t; mkdir -p output_bdd_$t
  cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_bdd_$t/; done
arm bdd_new ../cli/hyd_dd config_bdd_new.xml
arm bdd_on  ../cli/hyd_dd config_bdd_on.xml HYD_DEEP_DRAG=1
wait_arms
cmp_dirs    bdd_new bzg_new "HYD A  OFF-BRANCH (hyd_dd unset vs hyd_zg unset)"
want_differ bdd_on  bdd_new "HYD C  CONTROL =1 day vs unset -- MUST differ"
grep -h "HYD_DEEP_DRAG tau" bdd_on.log | head -1
touch DD_VERIFY_DONE
