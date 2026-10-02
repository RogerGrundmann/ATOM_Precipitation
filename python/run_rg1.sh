#!/bin/bash
# 2026-10-02: MC-RG probe -- why do the rain forests not convect with the mixed-layer parcel? Restart wb9 600 -> 640, wb9 knobs, cli/atm_rg.
# PRE-REGISTERED: Amazon/Congo parcel theta_e falls short of the env theta_es minimum by more than the ocean's (ocean ~0, forests
# -3..-6 K), and the shortfall is mostly the WARMER land free troposphere (T at 5 km ~2-3 K above the tropical ocean), not a
# weaker parcel (forest parcel theta_e within ~2 K of the ocean's).
set -u; cd "$(dirname "$0")"; rm -f RG1_DONE
mkdir output_rg1 || { touch RG1_DONE; exit 1; }
[ -e config_rg1.xml ] && { touch RG1_DONE; exit 1; }
ln -s ../output_wb9/atm_restart_0Ma_600.bin output_rg1/atm_restart_0Ma_600.bin
sed -e "s#output_wb9/#output_rg1/#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
    -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" -e "s#<nm>600<#<nm>640<#" config_wb9.xml > config_rg1.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_rg $(md5sum < ../cli/atm_rg | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_RH_OCEAN=0.82 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=1 \
    ../cli/atm_rg config_rg1.xml > rg1.log 2>&1
echo "rg1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' rg1.log)  $(date +%H:%M)"
grep -a "\[MC-RG\]" rg1.log | tail -6
touch RG1_DONE
