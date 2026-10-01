#!/bin/bash
# 2026-10-01: MC-LO probe -- the convection split by surface type ([MC-LO], ATM_MC_DIAG): why does wb7 rain 1000 mm/a over
# land and 226 over ocean (NASA 782 / 1056)? Restart 600 -> 640 from wb7's own checkpoint, current working branch
# (acb9c2c env), cli/atm_lo (-O2), 8 threads.
# PRE-REGISTERED (hypotheses, not predictions of numbers): ocean columns trigger less (lower active fraction) or have
# little CAPE / a dry sub-cloud layer; or they trigger but evaporate most of the rain (low P_conv / g_p).
set -u; cd "$(dirname "$0")"; rm -f LO1_DONE
mkdir output_lo1 || { touch LO1_DONE; exit 1; }
[ -e config_lo1.xml ] && { touch LO1_DONE; exit 1; }
ln -s ../output_wb7/atm_restart_0Ma_600.bin output_lo1/atm_restart_0Ma_600.bin
sed -e "s#output_wb7/#output_lo1/#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" -e "s#<nm>600<#<nm>640<#" config_wb7.xml > config_lo1.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_lo $(md5sum < ../cli/atm_lo | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ../cli/atm_lo config_lo1.xml > lo1.log 2>&1
echo "lo1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' lo1.log)  $(date +%H:%M)"
grep -a "\[MC-LO\]" lo1.log | tail -4
grep -a -i "land .*ocean" lo1.log | tail -1 | cut -c1-100
touch LO1_DONE
