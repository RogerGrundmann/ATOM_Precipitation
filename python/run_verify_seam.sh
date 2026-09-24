#!/bin/bash
# ATM_SEAM_PERIODIC (new 2026-09-24, default 0). SHORT CHECKS ONLY -- longer runs postponed by the user.
# 1 thread, nm = 20 from scratch, cli/atm_sp.
#   bsp_new  unset                                         -- MUST equal bwc_new (atm_wc unset), RUN_CONFIG aside
#   bsp_fo0  ATM_BC_SECOND_ORDER=0                         -- first-order reference
#   bsp_fo1  ATM_BC_SECOND_ORDER=0 ATM_SEAM_PERIODIC=1     -- MUST equal bsp_fo0: a no-op by arithmetic at first order
#   bsp_on   ATM_SEAM_PERIODIC=1                           -- CONTROL: MUST differ from bsp_new
# NOT tested here, and it cannot be at 20 iterations: whether it removes the k = 1 seam mode, which
# needs ~300 iterations to leave the noise (2026-09-21). That is the postponed 600-iteration arm.
set -u; cd "$(dirname "$0")"; rm -f SP_VERIFY_DONE
until [ -f WC_VERIFY_DONE ]; do sleep 30; done
. ./verify_lib.sh
for t in new fo0 fo1 on; do rm -rf output_bsp_$t; mkdir -p output_bsp_$t; done
arm bsp_new ../cli/atm_sp config_bsp_new.xml
arm bsp_fo0 ../cli/atm_sp config_bsp_fo0.xml ATM_BC_SECOND_ORDER=0
arm bsp_fo1 ../cli/atm_sp config_bsp_fo1.xml ATM_BC_SECOND_ORDER=0 ATM_SEAM_PERIODIC=1
arm bsp_on  ../cli/atm_sp config_bsp_on.xml  ATM_SEAM_PERIODIC=1
wait_arms
cmp_dirs    bsp_new bwc_new "ATM A  OFF-BRANCH (atm_sp unset vs atm_wc unset)"
cmp_dirs    bsp_fo1 bsp_fo0 "ATM B  NO-OP AT FIRST ORDER (=1 vs unset, both BC_SECOND_ORDER=0)"
want_differ bsp_on  bsp_new "ATM C  CONTROL =1 vs unset on the second-order default -- MUST differ"
touch SP_VERIFY_DONE
