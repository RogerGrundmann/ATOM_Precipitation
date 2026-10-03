#!/bin/bash
# 2026-10-03: wb19a / wb19b = wb14 + ATM_MC_T_ADD_LAND 2.0 / 3.0 K (user: land-only parcel excess). SCREENING: nm 220 from scratch, cli/atm_tl (-O2),
# 2 x 8 threads. One variable against wb14.
# wb14 at 200: 849 mm/a, r .458, sigma 3.74, land/ocean 194/1108, bands 2942/131/210/4.5; land 0-15 / 15-35 / 35-65 / 65-90 = 354 / 1 / 385 / 7;
# parcel theta_e - env theta_es min: Amazon -2.5 K, Congo -3.8, SE Asia -2.9, Sahara -14.4, Arabia -12.1 (ocean<15 +2.3). theta_e gains ~1.16 K per K.
# PRE-REGISTERED: margins rise by ~2.3 K (a) / ~3.5 K (b): a -- Amazon ~-0.2, Congo -1.5, SE Asia -0.6: weak forest rain, Amazon 0.5-3 mm/d, Congo 2-4
# (its 2.0 is stratiform), land 0-15 354 -> 500-1200; b -- Amazon +1, Congo -0.3, SE Asia +0.6: Amazon 3-10 mm/d, Congo 3-7, land 0-15 1200-3000.
# Both: Sahara, Arabia, Oman, W Australia, Namib stay 0 (margins <= -8 K); ocean unchanged (1108 +- 1 %); Horn stays dry (cap);
# land 15-35 < 100; r >= .44. RISK: hot moist land outside the listed regions (Sahel, India, N Australia) fires in single cells > 50 mm/d.
set -u; cd "$(dirname "$0")"; rm -f WB19_DONE
for t in wb19a wb19b; do
  mkdir output_$t || { touch WB19_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB19_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_tl $(md5sum < ../cli/atm_tl | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.805 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 ATM_MC_GATE_BLEND=1 ATM_RH_LAND_QCAP=1"
env OMP_NUM_THREADS=8 $K ATM_MC_T_ADD_LAND=2.0 ../cli/atm_tl config_wb19a.xml > wb19a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_MC_T_ADD_LAND=3.0 ../cli/atm_tl config_wb19b.xml > wb19b.log 2>&1 &
wait
for t in wb19a wb19b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean" $t.log | tail -1
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
  tail -1 output_$t/convergence.csv
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
touch WB19_DONE
