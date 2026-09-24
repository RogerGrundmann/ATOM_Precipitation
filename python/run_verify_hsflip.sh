#!/bin/bash
# BOTH-DIRECTIONS BYTE CHECK for ATM_HYDRO_SPLIT 0.0 -> 1.0 (default), 2026-09-24. 1 thread, nm = 20 from scratch.
#   bhs_old  cli/atm_prehs (= atm_wc2, HEAD c546072, split default 0), clean  = the old branch
#   bhs_off  cli/atm_hs1 (flipped), ATM_HYDRO_SPLIT=0                        -- MUST equal bhs_old
#   bhs_new  cli/atm_hs1, clean                                              = the new default
#   bhs_on   cli/atm_prehs, ATM_HYDRO_SPLIT=1.0                              -- MUST equal bhs_new
#   CONTROL  bhs_new vs bhs_off MUST differ. RUN_CONFIG.txt differs by banner token/path only.
set -u; cd "$(dirname "$0")"; rm -f HSFLIP_VERIFY_DONE
. ./verify_lib.sh
for t in old off new on; do rm -rf output_bhs_$t; mkdir -p output_bhs_$t; done
arm bhs_old ../cli/atm_prehs config_bhs_old.xml
arm bhs_off ../cli/atm_hs1   config_bhs_off.xml ATM_HYDRO_SPLIT=0
arm bhs_new ../cli/atm_hs1   config_bhs_new.xml
arm bhs_on  ../cli/atm_prehs config_bhs_on.xml  ATM_HYDRO_SPLIT=1.0
wait_arms
cmp_dirs    bhs_off bhs_old "A  RESTORE direction (new binary =0 vs pre-flip clean)"
cmp_dirs    bhs_on  bhs_new "B  FLIP direction    (pre-flip =1.0 vs new binary clean)"
want_differ bhs_new bhs_off "C  CONTROL new default vs off -- MUST differ"
for p in "off old" "on new"; do set -- $p
  diff <(sed "s#output_bhs_$1#X#g" output_bhs_$1/RUN_CONFIG.txt) <(sed "s#output_bhs_$2#X#g" output_bhs_$2/RUN_CONFIG.txt) | grep -o "HYDRO_SPLIT=[^ ]*" | tr '\n' ' '; echo; done
touch HSFLIP_VERIFY_DONE
