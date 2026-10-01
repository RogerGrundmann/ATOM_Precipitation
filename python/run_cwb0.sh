#!/bin/bash
# 2026-10-01: verification -- ColumnWaterBudget excludes ocean level 0 (the SST skin). Restart 600 -> 640 from
# output_rhst_115, current working branch (incl. ATM_Q_DIFF_FLUX via the env), cli/atm_cwb0 (-O2), 8 threads.
# PRE-REGISTERED: [CWB] evaporation row = the printed Evaporation_average (~706); NET falls by ~400; every other row
# within a few %.
set -u; cd "$(dirname "$0")"; rm -f CWB0_DONE
mkdir output_cwb0 || { touch CWB0_DONE; exit 1; }
[ -e config_cwb0.xml ] && { touch CWB0_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_cwb0/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_cwb0/#" \
    -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" \
    -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_cwb0.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_cwb0 $(md5sum < ../cli/atm_cwb0 | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_CWB_DIAG=1 \
    ../cli/atm_cwb0 config_cwb0.xml > cwb0.log 2>&1
echo "cwb0 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' cwb0.log)  $(date +%H:%M)"
touch CWB0_DONE
