#!/bin/bash
# BOTH-DIRECTIONS BYTE CHECK for ATM_MC_S_NDIM 0 -> 1 (default), 2026-09-24. 1 thread, nm = 20 from scratch.
set -u; cd "$(dirname "$0")"; rm -f MSF_VERIFY_DONE
. ./verify_lib.sh
for t in old off new on; do rm -rf output_msf_$t; mkdir -p output_msf_$t; done
arm msf_old ../cli/atm_pt482  config_msf_old.xml
arm msf_off ../cli/atm_msflip config_msf_off.xml ATM_MC_S_NDIM=0
arm msf_new ../cli/atm_msflip config_msf_new.xml
arm msf_on  ../cli/atm_pt482  config_msf_on.xml  ATM_MC_S_NDIM=1
wait_arms
cmp_dirs    msf_off msf_old "A  RESTORE direction"
cmp_dirs    msf_on  msf_new "B  FLIP direction"
want_differ msf_new msf_off "C  CONTROL -- MUST differ"
for p in "off old" "on new"; do set -- $p
  diff <(sed "s#output_msf_$1#X#g" output_msf_$1/RUN_CONFIG.txt) <(sed "s#output_msf_$2#X#g" output_msf_$2/RUN_CONFIG.txt) | grep -o "MC_S_NDIM=[^ ]*" | tr '\n' ' '; echo; done
touch MSF_VERIFY_DONE
