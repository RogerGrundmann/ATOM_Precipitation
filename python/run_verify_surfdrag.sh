#!/bin/bash
# OFF-BRANCH BYTE CHECK for ATM_SURF_DRAG_CONSISTENT (B.9, new 2026-09-24, default 0.0).
# 1 thread, nm = 20 from scratch.
#   bsd_old  cli/atm_rt (4305c587, HEAD aa47176 atmosphere), clean
#   bsd_new  cli/atm_sd (ee2ddaa4), unset  -- MUST equal bsd_old except RUN_CONFIG.txt (new banner token)
#   bsd_on   cli/atm_sd, =1.0              -- CONTROL: MUST differ
set -u; cd "$(dirname "$0")"; rm -f SD_VERIFY_DONE
. ./verify_lib.sh
for t in old new on; do rm -rf output_bsd_$t; mkdir -p output_bsd_$t; done
arm bsd_old ../cli/atm_rt config_bsd_old.xml
arm bsd_new ../cli/atm_sd config_bsd_new.xml
arm bsd_on  ../cli/atm_sd config_bsd_on.xml ATM_SURF_DRAG_CONSISTENT=1.0
wait_arms
cmp_dirs    bsd_new bsd_old "ATM A  OFF-BRANCH (atm_sd unset vs HEAD)"
diff <(sed 's#output_bsd_new#X#g' output_bsd_new/RUN_CONFIG.txt) <(sed 's#output_bsd_old#X#g' output_bsd_old/RUN_CONFIG.txt)
want_differ bsd_on  bsd_new "ATM C  CONTROL =1 vs unset -- MUST differ"
grep -h "surface drag coefficient" bsd_on.log | head -1
touch SD_VERIFY_DONE
