#!/bin/bash
# OFF-BRANCH BYTE CHECK for ATM_MC_QVD (default 0), 2026-09-25. 1 thread, nm = 20. old = cli/atm_bst (4272909), new = cli/atm_qvd.
set -u; cd "$(dirname "$0")"; rm -f QVD_VERIFY_DONE
. ./verify_lib.sh
for t in old new ctl; do mkdir output_qvd_$t || exit 1; done
arm qvd_old ../cli/atm_bst config_qvd_old.xml
arm qvd_new ../cli/atm_qvd config_qvd_new.xml
arm qvd_ctl ../cli/atm_qvd config_qvd_ctl.xml ATM_MC_QVD=2
wait_arms
cmp_dirs    qvd_new qvd_old "A  OFF BRANCH (new clean vs old)"
want_differ qvd_ctl qvd_new "C  CONTROL -- MUST differ"
diff <(sed "s#output_qvd_new#X#g" output_qvd_new/RUN_CONFIG.txt) <(sed "s#output_qvd_old#X#g" output_qvd_old/RUN_CONFIG.txt) | grep -o "MC_QVD=[^ ]*" | tr '\n' ' '; echo
touch QVD_VERIFY_DONE
