#!/bin/bash
# 2026-10-01: a panorama VTS for output_hsf1800b_on (that run had paraview_panorama_vts_flag=false). Restart from its
# atm_restart_0Ma_1800.bin with the SAME binary (cli/atm_hsf, -O0, 09-23) and the SAME env (ATM_RAD_TOPO=0
# ATM_UBUD_BALANCE=1 ATM_HYDRO_SPLIT=1.0), 2 iterations (1800 -> 1802, so the moist block refreshes the precipitation /
# cloud diagnostics), panorama at 1802. The .vts is the state 2 iterations (0.4 s) after the checkpoint.
set -u; cd "$(dirname "$0")"; rm -f VTS1800_DONE
mkdir output_vts1800 || { touch VTS1800_DONE; exit 1; }
[ -e config_vts1800.xml ] && { touch VTS1800_DONE; exit 1; }
ln -s ../output_hsf1800b_on/atm_restart_0Ma_1800.bin output_vts1800/atm_restart_0Ma_1800.bin
sed -e "s#output_hsf1800b_on/#output_vts1800/#" \
    -e "s#<checkpoint_save_iter>1800<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>1200<#<restart_from_iter>1800<#" \
    -e "s#<nm>1800<#<nm>1802<#" \
    -e "s#<paraview_panorama_vts_flag>false<#<paraview_panorama_vts_flag>true<#" \
    -e "s#<panorama_print>1200<#<panorama_print>1802<#" config_hsf1800b_on.xml > config_vts1800.xml
echo "start $(date +%H:%M)  atm_hsf $(md5sum < ../cli/atm_hsf | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_RAD_TOPO=0 ATM_UBUD_BALANCE=1 ATM_HYDRO_SPLIT=1.0 \
    ../cli/atm_hsf config_vts1800.xml > vts1800.log 2>&1
echo "vts1800 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' vts1800.log)  $(date +%H:%M)"
ls -la output_vts1800/*.vts 2>&1
touch VTS1800_DONE
