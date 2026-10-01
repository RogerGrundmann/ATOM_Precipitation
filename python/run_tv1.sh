#!/bin/bash
# 2026-10-01: MC-TV probe -- coeff_MC_t / coeff_MC_vel (L_atm, 2.5 % of L/u_0): convective energy vs L*P_conv, momentum
# tendencies in m/s/day at the scheme value and as applied, cap census, flux ingredients. Restart 600 -> 640 from
# output_rhst_115, current working branch (all four water fixes), cli/atm_tv (-O2), 8 threads.
# PRE-REGISTERED: applied heating ~2.5 % of L*P_conv (~30 W/m2); scheme |MC_w| p50 >> 1-3 m/s/day observed CMT.
set -u; cd "$(dirname "$0")"; rm -f TV1_DONE
mkdir output_tv1 || { touch TV1_DONE; exit 1; }
[ -e config_tv1.xml ] && { touch TV1_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_tv1/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_tv1/#" \
    -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" \
    -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_tv1.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_tv1 $(md5sum < ../cli/atm_tv | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 \
    ../cli/atm_tv config_tv1.xml > tv1.log 2>&1
echo "tv1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' tv1.log)  $(date +%H:%M)"
touch TV1_DONE
