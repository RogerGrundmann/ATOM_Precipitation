#!/bin/bash
# OFF-BRANCH BYTE CHECK for ATM_MC_BASE_SAT=2 and ATM_MC_GP_AREA (both default 0), 2026-09-26. 1 thread, nm = 20.
# old = cli/atm_dth, new = cli/atm_gpa. Two controls on the alf_tie stack, one per change. Queued behind qth_on.
set -u; cd "$(dirname "$0")"; rm -f DGPA_VERIFY_DONE
until [ -f QTH100_DONE ]; do sleep 60; done
. ./verify_lib.sh
for t in old new stk ctlb ctlc; do mkdir output_dgpa_$t || { touch DGPA_VERIFY_DONE; exit 1; }; done
STK="ATM_MC_SGZ=1 ATM_MC_ENTR=1.0e-4 ATM_MC_QVD=2 ATM_MC_ALF1=5.44e-4"
arm dgpa_old  ../cli/atm_dth config_dgpa_old.xml
arm dgpa_new  ../cli/atm_gpa config_dgpa_new.xml
arm dgpa_stk  ../cli/atm_gpa config_dgpa_stk.xml  $STK ATM_MC_BASE_SAT=1
arm dgpa_ctlb ../cli/atm_gpa config_dgpa_ctlb.xml $STK ATM_MC_BASE_SAT=2
arm dgpa_ctlc ../cli/atm_gpa config_dgpa_ctlc.xml $STK ATM_MC_BASE_SAT=1 ATM_MC_GP_AREA=1
wait_arms
cmp_dirs    dgpa_new  dgpa_old "A  OFF BRANCH (new clean vs old)"
want_differ dgpa_ctlb dgpa_stk "C  CONTROL (b) BASE_SAT=2 -- MUST differ"
want_differ dgpa_ctlc dgpa_stk "C  CONTROL (c) GP_AREA=1 -- MUST differ"
diff <(sed "s#output_dgpa_new#X#g" output_dgpa_new/RUN_CONFIG.txt) <(sed "s#output_dgpa_old#X#g" output_dgpa_old/RUN_CONFIG.txt) | grep -o "MC_GP_AREA=[^ ]*" | tr '\n' ' '; echo
touch DGPA_VERIFY_DONE
