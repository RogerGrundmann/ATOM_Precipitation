#!/bin/bash
# BOTH-DIRECTIONS BYTE CHECK, 2026-09-23: c43 = 4/3, c13 = 1/3 (second-order Neumann) made the
# DEFAULT in both models, at the user's instruction. 1 thread, nm = 20 from scratch.
#   atmosphere: cli/atm_hs3 (md5 5d0a3b11, same source with the first-order default) vs cli/atm_bc2 (514e9c17)
#   ocean:      cli/hyd_prebc (34a65bd8) vs cli/hyd_bc2 (63d3c910)
# A = RESTORE (new binary =0 vs old clean), B = FLIP (old =1 vs new clean), C = CONTROL must differ.
set -u; cd "$(dirname "$0")"; rm -f BC2_VERIFY_DONE
. ./verify_lib.sh
for t in old off new on; do rm -rf output_bca_$t output_bch_$t; mkdir -p output_bca_$t output_bch_$t; done
arm bca_old ../cli/atm_hs3 config_bca_old.xml
arm bca_off ../cli/atm_bc2 config_bca_off.xml ATM_BC_SECOND_ORDER=0
arm bca_new ../cli/atm_bc2 config_bca_new.xml
arm bca_on  ../cli/atm_hs3 config_bca_on.xml  ATM_BC_SECOND_ORDER=1
arm bch_old ../cli/hyd_prebc config_bch_old.xml
arm bch_off ../cli/hyd_bc2   config_bch_off.xml HYD_BC_SECOND_ORDER=0
arm bch_new ../cli/hyd_bc2   config_bch_new.xml
arm bch_on  ../cli/hyd_prebc config_bch_on.xml  HYD_BC_SECOND_ORDER=1
wait_arms
cmp_dirs    bca_off bca_old "ATM A  RESTORE"
cmp_dirs    bca_on  bca_new "ATM B  FLIP"
want_differ bca_new bca_off "ATM C  CONTROL"
cmp_dirs    bch_off bch_old "HYD A  RESTORE"
cmp_dirs    bch_on  bch_new "HYD B  FLIP"
want_differ bch_new bch_off "HYD C  CONTROL"
touch BC2_VERIFY_DONE
