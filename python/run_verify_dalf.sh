#!/bin/bash
# OFF-BRANCH BYTE CHECK for ATM_MC_ALF1 (default 0.05 = shipped), 2026-09-26. 1 thread, nm = 20 from scratch.
# old = cli/atm_dqh, new = cli/atm_alf. The CONTROL is run on the qvd2 stack (SGZ+ENTR 1e-4+BASE_SAT+QVD=2),
# because on the shipped branch e_p is ~0 and the knob could not show; stk = that stack, ctl = stack + ALF1.
set -u; cd "$(dirname "$0")"; rm -f DALF_VERIFY_DONE
. ./verify_lib.sh
for t in old new stk ctl; do mkdir output_dalf_$t || exit 1; done
STK="ATM_MC_SGZ=1 ATM_MC_ENTR=1.0e-4 ATM_MC_BASE_SAT=1 ATM_MC_QVD=2"
arm dalf_old ../cli/atm_dqh config_dalf_old.xml
arm dalf_new ../cli/atm_alf config_dalf_new.xml
arm dalf_stk ../cli/atm_alf config_dalf_stk.xml $STK
arm dalf_ctl ../cli/atm_alf config_dalf_ctl.xml $STK ATM_MC_ALF1=5.44e-4
wait_arms
cmp_dirs    dalf_new dalf_old "A  OFF BRANCH (new clean vs old)"
want_differ dalf_ctl dalf_stk "C  CONTROL -- MUST differ"
diff <(sed "s#output_dalf_new#X#g" output_dalf_new/RUN_CONFIG.txt) <(sed "s#output_dalf_old#X#g" output_dalf_old/RUN_CONFIG.txt) | grep -o "MC_ALF1=[^ ]*" | tr '\n' ' '; echo
touch DALF_VERIFY_DONE
