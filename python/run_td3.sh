#!/bin/bash
# 2026-10-01: DIFF-LEAK discriminator -- the radial 2nd-derivative diffusion row (+222 mm/a in td2) re-summed with
# rho * centred thickness, const rho * forward thickness, const rho * centred thickness. Restart 600 -> 640 from
# output_rhst_115, current working branch, cli/atm_td3 (-O2), 8 threads.
# PRE-REGISTERED: the centred-thickness sums stay large (the operator does not telescope on the stretched grid);
# the const-rho sums change sign or size only modestly.
set -u; cd "$(dirname "$0")"; rm -f TD3_DONE
mkdir output_td3 || { touch TD3_DONE; exit 1; }
[ -e config_td3.xml ] && { touch TD3_DONE; exit 1; }
ln -s ../output_rhst_115/atm_restart_0Ma_600.bin output_td3/atm_restart_0Ma_600.bin
sed -e "s#output_rhst_115/#output_td3/#" \
    -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" \
    -e "s#<nm>600<#<nm>640<#" config_rhst_115.xml > config_td3.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_td3 $(md5sum < ../cli/atm_td3 | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_CWB_DIAG=1 \
    ../cli/atm_td3 config_td3.xml > td3.log 2>&1
echo "td3 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' td3.log)  $(date +%H:%M)"
touch TD3_DONE
