#!/bin/bash
# 2026-10-02: ATM_TEQ_WTG=1 probe -- restart 600 -> 640 from wb7, working branch, cli/atm_wtg (-O2), 8 threads.
# Control: lo2 (same restart, knob absent). t_eq is rebuilt from the WTG baseline at iteration 620 (teq_refresh_stride).
# PRE-REGISTERED: |lat|<30 land-ocean T and t_eq at 5/8 km fall from +2.3..2.8 K to < 0.5 K; ocean CAPE rises toward
# land (638 -> >900 J/kg), land CAPE falls (1252 -> <1000); ocean P_conv rises, land P_conv falls.
set -u; cd "$(dirname "$0")"; rm -f LO3_DONE
mkdir output_lo3 || { touch LO3_DONE; exit 1; }
[ -e config_lo3.xml ] && { touch LO3_DONE; exit 1; }
ln -s ../output_wb7/atm_restart_0Ma_600.bin output_lo3/atm_restart_0Ma_600.bin
sed -e "s#output_wb7/#output_lo3/#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" -e "s#<nm>600<#<nm>640<#" config_wb7.xml > config_lo3.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_wtg $(md5sum < ../cli/atm_wtg | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_TEQ_WTG=1 ../cli/atm_wtg config_lo3.xml > lo3.log 2>&1
echo "lo3 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' lo3.log)  $(date +%H:%M)"
grep -a "TEQ-WTG" lo3.log
grep -a "\[MC-LO\]" lo3.log | tail -6
grep -a -i "land .*ocean" lo3.log | tail -1 | cut -c1-100
touch LO3_DONE
