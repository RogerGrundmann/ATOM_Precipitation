#!/bin/bash
# 2026-10-04: wb26a / wb26b = wb25b + ATM_RH_LAND_EAST 0.5 / 1.0 (user). wb25b = working branch + ML parcel stack + ATM_RH_LAND=2 + QCAP +
# ATM_MC_GATE_BLEND=3 + ATM_MC_GP_AREA=1 + ATM_RH_OCEAN=0.80 + ATM_MC_T_ADD_LAND=2.3. ONE VARIABLE against wb25b.
# SCREENING: nm 220 from scratch, cli/atm_g3 (-O2, = cli/atm at 759f496), 2 x 8 threads.
# wb25b: 912 mm/a, r .562, sigma 2.85, bands 3001/290/204/4.5, land/ocean 443/1098, land 0-15 / 15-35 / 35-65 / 65-90 = 1641 / 2 / 380 / 7,
# Amazon 8.25, Congo 3.12, India 0.41 mm/d; S China, SE US, E Australia, SE Africa, S Brazil 0 (NASA 3.87 / 3.89 / 1.86 / 2.22 / 4.55 mm/d).
# Earlier: wb17 (EAST 0.5 on wb14, no land excess, blend 1) was a NULL; wb16 (EAST 1.0 + SIGMA_SFC) land 15-35 853, SE Africa 4.64, S China 1.91 mm/d.
# PRE-REGISTERED: with the 2.3 K land excess and the undiluted gate the moistened east sides now convect.
# a (0.5): land 15-35 2 -> 20-150 mm/a, at most two of the five east-side regions above 0.3 mm/d; b (1.0): land 15-35 150-500 (NASA 643),
# S China / SE Africa / S Brazil / E Australia 0.5-4 mm/d, SE US < 0.5 (35N, outside the taper). Both: Sahara, Arabia, Oman, Namib,
# W Australia 0; ocean 1098 +-1 %; land 0-15 1641 +-5 %, Amazon / Congo within 5 %; 35-65, 65-90 unchanged; no land cell > 40 mm/d;
# r .555-.58 (b highest); global 912 -> a 915-935, b 930-1000. RISK: a null again at 0.5; QCAP caps the added land vapour (b far below range).
set -u; cd "$(dirname "$0")"; rm -f WB26_DONE
for t in wb26a wb26b; do
  mkdir output_$t || { touch WB26_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB26_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_g3 $(md5sum < ../cli/atm_g3 | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 ATM_MC_GATE_BLEND=3 ATM_MC_GP_AREA=1 ATM_RH_OCEAN=0.80 ATM_MC_T_ADD_LAND=2.3 ATM_RH_LAND_QCAP=1"
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=0.5 ../cli/atm_g3 config_wb26a.xml > wb26a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.0 ../cli/atm_g3 config_wb26b.xml > wb26b.log 2>&1 &
wait
for t in wb26a wb26b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "RUN CONFIG" $t.log | grep -o "MC_GP_AREA=[^ ]*\|MC_GATE_BLEND=[^ ]*\|RH_OCEAN=[^ ]*\|MC_T_ADD_LAND=[^ ]*\|RH_LAND_EAST=[^ ]*" | sort -u | tr '\n' ' '; echo
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean\|P_rain mean" $t.log | tail -2
  grep -a "\[MC-GATE\]" $t.log | tail -2 | cut -c1-260
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
  tail -1 output_$t/convergence.csv
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
touch WB26_DONE
