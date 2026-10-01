#!/bin/bash
# 2026-10-01: why is c_u = 0? [MC-Q] updraftRecurrence census (|M_u| bins vs coeff_recurr = 0.1, parcel
# supersaturation where the condensation step runs / is skipped, c_u written). Restart 600 -> 640 from output_rhst_115,
# working branch, cli/atm_cuz (-O2), 8 threads. PRE-REGISTERED: |M_u| <= 0.1 in nearly every cell, so the step never runs.
set -u; cd "$(dirname "$0")"; rm -f CUZ_DONE
mkdir output_cuz || { touch CUZ_DONE; exit 1; }
[ -e config_cuz.xml ] && { touch CUZ_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_cuz/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_cuz/#" \
    -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" \
    -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_cuz.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_cuz $(md5sum < ../cli/atm_cuz | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 \
    ../cli/atm_cuz config_cuz.xml > cuz.log 2>&1
echo "cuz exit $?  NaN $(grep -c 'NaN/Inf DETECTED' cuz.log)  $(date +%H:%M)"
touch CUZ_DONE
