#!/bin/bash
# 2026-10-01: MC-TV pair probe -- ATM_MC_T_COEFF=1 + ATM_MC_UV_DETRAIN=1 + ATM_MC_VEL_COEFF=1. Restart 600 -> 640 from
# output_rhst_115, working branch, cli/atm_tv3 (-O2), 8 threads. Compare tv2 (VEL_COEFF off) and tv1 (all off).
# PRE-REGISTERED: applied |MC_w| = the scheme value (p50/p90/p99 ~0.1/6/24 m/s/day); exit 0, no NaN, max winds within
# a few % of tv2; w-budget drag_conv rises ~40x from tv2 but stays below Coriolis in the median.
set -u; cd "$(dirname "$0")"; rm -f TV3_DONE
mkdir output_tv3 || { touch TV3_DONE; exit 1; }
[ -e config_tv3.xml ] && { touch TV3_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_tv3/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_tv3/#" \
    -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" \
    -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_tv3.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_tv3 $(md5sum < ../cli/atm_tv3 | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_MC_T_COEFF=1 ATM_MC_UV_DETRAIN=1 ATM_MC_VEL_COEFF=1 \
    ../cli/atm_tv3 config_tv3.xml > tv3.log 2>&1
echo "tv3 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' tv3.log)  $(date +%H:%M)"
touch TV3_DONE
