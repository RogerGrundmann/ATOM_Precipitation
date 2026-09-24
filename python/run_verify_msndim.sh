#!/bin/bash
# OFF-BRANCH BYTE CHECK for ATM_MC_S_NDIM (B.10 remainder, new 2026-09-24, default 0). 1 thread, nm = 20.
#   bms_old  cli/atm_bcd (tree just before this knob), clean
#   bms_new  cli/atm_msn unset -- MUST equal bms_old except RUN_CONFIG.txt (banner token)
#   bms_on   cli/atm_msn =1    -- CONTROL: MUST differ
set -u; cd "$(dirname "$0")"; rm -f MSN_VERIFY_DONE
. ./verify_lib.sh
for t in old new on; do rm -rf output_bms_$t; mkdir -p output_bms_$t; done
arm bms_old ../cli/atm_bcd config_bms_old.xml
arm bms_new ../cli/atm_msn config_bms_new.xml
arm bms_on  ../cli/atm_msn config_bms_on.xml ATM_MC_S_NDIM=1
wait_arms
cmp_dirs    bms_new bms_old "ATM A  OFF-BRANCH (atm_msn unset vs atm_bcd)"
diff <(sed 's#output_bms_new#X#g' output_bms_new/RUN_CONFIG.txt) <(sed 's#output_bms_old#X#g' output_bms_old/RUN_CONFIG.txt) | grep -o "MC_S_NDIM=[^ ]*"
want_differ bms_on  bms_new "ATM C  CONTROL =1 vs unset -- MUST differ"
touch MSN_VERIFY_DONE
