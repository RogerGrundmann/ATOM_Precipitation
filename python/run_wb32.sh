#!/bin/bash
# 2026-10-04: wb32a / wb32b = the wb31a stack + ATM_RH_OCEAN_ML 1000 / 1500 m (user: close STORM first). wb31a stack = working branch +
# ATM_RH_OCEAN=0.82 + ATM_MC_MB_SAT_OCEAN = ATM_MC_MB_SAT_LAND = 0.045 (tropics within 10 % of NASA). ONE VARIABLE against wb31a.
# SCREENING: nm 220 from scratch, cli/atm_oml (-O2, 211fd08 + the knob), 2 x 8 threads; the -O0 byte check (run_vmoml.sh) runs alongside.
# wb31a: 750 mm/a (-23 %), r .514, sigma 1.40, bands 2036/509/234/4.5; ocean 35-65 171 (NASA 1107), land 35-65 401 (643), ocean 15-35 709, ocean 65-90 2.
# Diagnosis (wb27, 87E): rain forms at 1.5-4 km and evaporates in the cloud-free layer below ~1.2 km -- 30 % reaches the sea at 40-56S, 6 % at 36S.
# The knob: constant q from the surface up to <depth> over ocean poleward of 30 deg (full from 40), RH capped at 0.98.
# STORM pass criteria (proposed): ocean 35-65 > 500, land 35-65 <= 770, r >= .55.
# PRE-REGISTERED: a (1000 m) -- ocean 35-65 171 -> 350-600, 35-65 band 234 -> 360-540, global 750 -> 800-860; b (1500 m, layer joins the deck) --
# ocean 35-65 450-800, band 430-680, global 830-910. Both: land 35-65 401 +-3 %, tropics (0-15, land 0-15, Amazon, Congo, ocean<=30 distribution)
# unchanged within 1 %, ocean 15-35 709 -> 720-800 (the 30-40 deg taper), ocean 65-90 2 -> 5-40, r rises to .53-.57, sigma 1.3-1.4, no NaN.
# RISK: a boundary layer at RH 0.98 against H_crit 0.7-0.9 is nearly overcast -- drizzle could overshoot (ocean 35-65 > 1500).
set -u; cd "$(dirname "$0")"; rm -f WB32_DONE
for t in wb32a wb32b; do
  mkdir output_$t || { touch WB32_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB32_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_oml $(md5sum < ../cli/atm_oml | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.82 ATM_MC_MB_SAT_OCEAN=0.045 ATM_MC_MB_SAT_LAND=0.045"
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN_ML=1000 ../cli/atm_oml config_wb32a.xml > wb32a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN_ML=1500 ../cli/atm_oml config_wb32b.xml > wb32b.log 2>&1 &
wait
for t in wb32a wb32b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  echo "banner diff vs wb31a (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb31a.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean\|P_rain mean" $t.log | tail -2
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-RG\] ocean<15" $t.log | tail -1 | cut -c1-260
  tail -1 output_$t/convergence.csv
  python3 oceanb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk | grep -E "India|Sahara|Arabia|Amazon|Congo|global max"
  python3 tropdist.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk | head -1
  python3 storm87.py output_$t/0Ma_smooth_Atm_zonal_87_220.vtk 2>/dev/null | grep -A8 "OCEAN 40-56S" | cut -c1-150
done
touch WB32_DONE
