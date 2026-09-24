#!/bin/bash
# BOTH-DIRECTIONS BYTE CHECK for ATM_SEAM_PERIODIC 0 -> 1 (default), 2026-09-24. 1 thread, nm = 20 from scratch.
#   bsf_old  cli/atm_fx (pre-flip, seam default 0), clean = old branch
#   bsf_off  cli/atm_spflip, ATM_SEAM_PERIODIC=0          -- MUST equal bsf_old
#   bsf_new  cli/atm_spflip, clean                        = new default
#   bsf_on   cli/atm_fx, ATM_SEAM_PERIODIC=1              -- MUST equal bsf_new
#   CONTROL  bsf_new vs bsf_off MUST differ.
set -u; cd "$(dirname "$0")"; rm -f SPFLIP_VERIFY_DONE
. ./verify_lib.sh
for t in old off new on; do rm -rf output_bsf_$t; mkdir -p output_bsf_$t; done
arm bsf_old ../cli/atm_fx     config_bsf_old.xml
arm bsf_off ../cli/atm_spflip config_bsf_off.xml ATM_SEAM_PERIODIC=0
arm bsf_new ../cli/atm_spflip config_bsf_new.xml
arm bsf_on  ../cli/atm_fx     config_bsf_on.xml  ATM_SEAM_PERIODIC=1
wait_arms
cmp_dirs    bsf_off bsf_old "A  RESTORE direction"
cmp_dirs    bsf_on  bsf_new "B  FLIP direction"
want_differ bsf_new bsf_off "C  CONTROL -- MUST differ"
for p in "off old" "on new"; do set -- $p
  diff <(sed "s#output_bsf_$1#X#g" output_bsf_$1/RUN_CONFIG.txt) <(sed "s#output_bsf_$2#X#g" output_bsf_$2/RUN_CONFIG.txt) | grep -o "SEAM_PERIODIC=[^ ]*" | tr '\n' ' '; echo; done
touch SPFLIP_VERIFY_DONE
