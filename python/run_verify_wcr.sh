#!/bin/bash
# BOTH-DIRECTIONS BYTE CHECK: ATM_WATER_CLOSURE default 1 -> 0 (2026-09-25). 1 thread, nm = 20 from scratch.
# old = cli/atm_sgz (23bee64, closure ON by default), new = cli/atm_wcr (closure OFF by default).
set -u; cd "$(dirname "$0")"; rm -f WCR_VERIFY_DONE
. ./verify_lib.sh
for t in old off new on; do mkdir output_wcr_$t || exit 1; done
arm wcr_old ../cli/atm_sgz config_wcr_old.xml ATM_WATER_CLOSURE=0
arm wcr_off ../cli/atm_wcr config_wcr_off.xml
arm wcr_new ../cli/atm_wcr config_wcr_new.xml ATM_WATER_CLOSURE=1
arm wcr_on  ../cli/atm_sgz config_wcr_on.xml
wait_arms
cmp_dirs    wcr_off wcr_old "A  RESTORE: new clean vs old ATM_WATER_CLOSURE=0 -- MUST MATCH"
cmp_dirs    wcr_new wcr_on  "B  FLIP: new ATM_WATER_CLOSURE=1 vs old clean  -- MUST MATCH"
want_differ wcr_off wcr_on "C  CONTROL new clean vs old clean -- MUST differ"
for p in "off old" "new on"; do set -- $p
  diff <(sed "s#output_wcr_$1#X#g" output_wcr_$1/RUN_CONFIG.txt) <(sed "s#output_wcr_$2#X#g" output_wcr_$2/RUN_CONFIG.txt) | grep -o "WATER_CLOSURE=[^ ]*" | tr '\n' ' '; echo; done
touch WCR_VERIFY_DONE
