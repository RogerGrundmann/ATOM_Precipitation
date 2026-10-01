#!/bin/bash
# 2026-10-01: MC-LO t_eq probe -- T and the prescribed relaxation target t_eq at 2/5/8/11 km, land minus ocean (global, |lat|<30).
# Restart 600 -> 640 from wb7, working branch, cli/atm_lo2 (-O2), 8 threads.
# PRE-REGISTERED: t_eq carries the same ~3 K land-ocean contrast aloft that T has (the relaxation holds T to it).
set -u; cd "$(dirname "$0")"; rm -f LO2_DONE
mkdir output_lo2 || { touch LO2_DONE; exit 1; }
[ -e config_lo2.xml ] && { touch LO2_DONE; exit 1; }
ln -s ../output_wb7/atm_restart_0Ma_600.bin output_lo2/atm_restart_0Ma_600.bin
sed -e "s#output_wb7/#output_lo2/#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" -e "s#<nm>600<#<nm>640<#" config_wb7.xml > config_lo2.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_lo $(md5sum < ../cli/atm_lo2 | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ../cli/atm_lo2 config_lo2.xml > lo2.log 2>&1
echo "lo2 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' lo2.log)  $(date +%H:%M)"
grep -a "\[MC-LO\]" lo2.log | tail -6
grep -a -i "land .*ocean" lo2.log | tail -1 | cut -c1-100
touch LO2_DONE
