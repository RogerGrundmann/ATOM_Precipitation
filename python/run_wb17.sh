#!/bin/bash
# 2026-10-03: wb17 = wb14 + ATM_RH_LAND_EAST=0.5, ATM_RH_SIGMA_SFC off (user). SCREENING arm: nm 220 from scratch (user 10-03: runs on this stack
# are within ~3 % of their iteration-600 values at 200), cli/atm_ea (-O2), 8 threads. One variable against wb14.
# wb14 at 200: 849 mm/a, r .458, sigma 3.74, land/ocean 194/1108, bands 2942/131/210/4.5; land 0-15 / 15-35 / 35-65 / 65-90 = 354 / 1 / 385 / 7.
# wb16 (EAST 1.0 + SIGMA_SFC, at 200): land 15-35 853 (wb15 without EAST: 198), SE Africa 4.64 mm/d (NASA 2.22), S China 1.91 (3.87).
# PRE-REGISTERED: land 15-35 1 -> 100-350 mm/a (NASA 643), from S China / SE Africa / S Brazil / E Australia; India and SE US still < 0.5 mm/d;
# Sahara, Arabia, Oman, Namib, W Australia stay 0; land 0-15, 35-65, 65-90 and the ocean unchanged against wb14 (within 3 %); global land 194 -> 230-320;
# r within +-0.02 of wb14's; no land cell above 40 mm/d.
set -u; cd "$(dirname "$0")"; rm -f WB17_DONE
mkdir output_wb17 || { touch WB17_DONE; exit 1; }
[ -e config_wb17.xml ] && { touch WB17_DONE; exit 1; }
sed -e "s#output_wb7/#output_wb17/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_wb17.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_ea $(md5sum < ../cli/atm_ea | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.805 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 \
    ATM_MC_GATE_BLEND=1 ATM_RH_LAND_QCAP=1 ATM_RH_LAND_EAST=0.5 ../cli/atm_ea config_wb17.xml > wb17.log 2>&1
t=wb17
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
grep -a "RH-LAND" $t.log | head -3
grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
grep -a "P_conv mean" $t.log | tail -1
grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
tail -1 output_$t/convergence.csv
python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
touch WB17_DONE
