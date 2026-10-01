#!/bin/bash
# 2026-10-01: MC-TV knob probe -- ATM_MC_T_COEFF=1 + ATM_MC_UV_DETRAIN=1. Restart 600 -> 640 from output_rhst_115,
# working branch, cli/atm_tv2 (-O2), 8 threads. Compare tv1 (both off, same window).
# PRE-REGISTERED: heating as applied = the scheme value ~ L*P_conv (~25 W/m2); |w_u - w| p50 falls from 15.5 to a few
# m/s, the scheme |MC_w| p50 from 112 toward 1-10 m/s/day, the cap census toward 0.
set -u; cd "$(dirname "$0")"; rm -f TV2_DONE
mkdir output_tv2 || { touch TV2_DONE; exit 1; }
[ -e config_tv2.xml ] && { touch TV2_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_tv2/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_tv2/#" \
    -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" \
    -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_tv2.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_tv2 $(md5sum < ../cli/atm_tv2 | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_MC_T_COEFF=1 ATM_MC_UV_DETRAIN=1 \
    ../cli/atm_tv2 config_tv2.xml > tv2.log 2>&1
echo "tv2 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' tv2.log)  $(date +%H:%M)"
touch TV2_DONE
