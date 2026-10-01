#!/bin/bash
# 2026-10-01: MC-Q-LEAK probe -- [MC-Q] rows split the convective moisture tendency (flux div, c_u, e_d, e_l, e_p, cap)
# against -P_conv. Restart 600 -> 640 from output_rhst_115 (working branch), cli/atm_mcq (-O2), 8 threads.
# PRE-REGISTERED: applied column MC_q ~ +383 x (1/0.713) physical; the leak sits in c_u - e_l - g_p (condensate
# neither rained nor re-evaporated) and/or the MCq_max cap; flux div ~0.
set -u; cd "$(dirname "$0")"; rm -f MCQ_DONE
mkdir output_mcq || { touch MCQ_DONE; exit 1; }
[ -e config_mcq.xml ] && { touch MCQ_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_mcq/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_mcq/#" \
    -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" \
    -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_mcq.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_mcq $(md5sum < ../cli/atm_mcq | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_SNOW_DIAG=1 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_SR_DIAG=1 \
    ../cli/atm_mcq config_mcq.xml > mcq.log 2>&1
echo "mcq exit $?  NaN $(grep -c 'NaN/Inf DETECTED' mcq.log)  $(date +%H:%M)"
touch MCQ_DONE
