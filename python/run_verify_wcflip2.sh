#!/bin/bash
# BOTH-DIRECTIONS BYTE CHECK: ATM_WATER_CLOSURE 0 -> 1 default, WITH ATM_SATADJ_FADE=2 folded in
# (2026-09-24). Replaces run_verify_wcflip.sh (stopped before it ran). 1 thread, nm = 20 from scratch.
#   bw2_old  cli/atm_prewc (closure default 0, fade default 0), clean
#   bw2_off  cli/atm_wcf2, ATM_WATER_CLOSURE=0                          -- MUST equal bw2_old
#   bw2_new  cli/atm_wcf2, clean
#   bw2_on   cli/atm_prewc, ATM_WATER_CLOSURE=1 ATM_SATADJ_FADE=2       -- MUST equal bw2_new
set -u; cd "$(dirname "$0")"; rm -f WCFLIP2_VERIFY_DONE
until [ -f SPFLIP_VERIFY_DONE ] && [ -f FXD3_DONE ]; do sleep 30; done
. ./verify_lib.sh
for t in old off new on; do rm -rf output_bw2_$t; mkdir -p output_bw2_$t; done
arm bw2_old ../cli/atm_prewc config_bw2_old.xml
arm bw2_off ../cli/atm_wcf2  config_bw2_off.xml ATM_WATER_CLOSURE=0
arm bw2_new ../cli/atm_wcf2  config_bw2_new.xml
arm bw2_on  ../cli/atm_prewc config_bw2_on.xml  ATM_WATER_CLOSURE=1 ATM_SATADJ_FADE=2
wait_arms
cmp_dirs    bw2_off bw2_old "A  RESTORE direction"
cmp_dirs    bw2_on  bw2_new "B  FLIP direction (closure + fade 2)"
want_differ bw2_new bw2_off "C  CONTROL -- MUST differ"
for p in "off old" "on new"; do set -- $p
  diff <(sed "s#output_bw2_$1#X#g" output_bw2_$1/RUN_CONFIG.txt) <(sed "s#output_bw2_$2#X#g" output_bw2_$2/RUN_CONFIG.txt) | grep -o "WATER_CLOSURE=[^ ]*\|SATADJ_FADE=[^ ]*" | tr '\n' ' '; echo; done
touch WCFLIP2_VERIFY_DONE
