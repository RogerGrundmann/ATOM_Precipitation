#!/bin/bash
# BOTH-DIRECTIONS BYTE CHECK for ATM_RH_MIN_PTOP 490 -> 482 (default), 2026-09-24. 1 thread, nm = 20.
set -u; cd "$(dirname "$0")"; rm -f P82_VERIFY_DONE
. ./verify_lib.sh
for t in old off new on; do rm -rf output_p82_$t; mkdir -p output_p82_$t; done
arm p82_old ../cli/atm_prept482 config_p82_old.xml
arm p82_off ../cli/atm_pt482    config_p82_off.xml ATM_RH_MIN_PTOP=490
arm p82_new ../cli/atm_pt482    config_p82_new.xml
arm p82_on  ../cli/atm_prept482 config_p82_on.xml  ATM_RH_MIN_PTOP=482
wait_arms
cmp_dirs    p82_off p82_old "A  RESTORE direction (=490 vs pre-flip)"
cmp_dirs    p82_on  p82_new "B  FLIP direction (pre-flip =482 vs new)"
want_differ p82_new p82_off "C  CONTROL -- MUST differ"
for p in "off old" "on new"; do set -- $p
  diff <(sed "s#output_p82_$1#X#g" output_p82_$1/RUN_CONFIG.txt) <(sed "s#output_p82_$2#X#g" output_p82_$2/RUN_CONFIG.txt) | grep -o "RH_MIN_PTOP=[^ ]*" | tr '\n' ' '; echo; done
touch P82_VERIFY_DONE
