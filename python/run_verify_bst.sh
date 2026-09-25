#!/bin/bash
# OFF-BRANCH BYTE CHECK for ATM_MC_BASE_SAT (default 0), 2026-09-25. 1 thread, nm = 20 from scratch.
# old = cli/atm_dtv (7795121), new = cli/atm_bst. Control: the knob WITH SGZ + ENTR (BASE_SAT alone may be
# near-inert where the parcel is already saturated at base).
set -u; cd "$(dirname "$0")"; rm -f BST_VERIFY_DONE
. ./verify_lib.sh
for t in old new ctl; do mkdir output_bst_$t || exit 1; done
arm bst_old ../cli/atm_dtv config_bst_old.xml
arm bst_new ../cli/atm_bst config_bst_new.xml
arm bst_ctl ../cli/atm_bst config_bst_ctl.xml ATM_MC_BASE_SAT=1
wait_arms
cmp_dirs    bst_new bst_old "A  OFF BRANCH (new clean vs old)"
want_differ bst_ctl bst_new "C  CONTROL -- MUST differ"
diff <(sed "s#output_bst_new#X#g" output_bst_new/RUN_CONFIG.txt) <(sed "s#output_bst_old#X#g" output_bst_old/RUN_CONFIG.txt) | grep -o "MC_BASE_SAT=[^ ]*" | tr '\n' ' '; echo
touch BST_VERIFY_DONE
