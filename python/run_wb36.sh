#!/bin/bash
# 2026-10-05: wb36a / wb36b = wb35b base + ATM_RH_LAND_EAST 1.0 / 1.4 (user approved 10-05; LAND-POLAR step 2, subtropical land).
# SCREENING: nm 220 from scratch, cli/atm_sgl (-O2, md5 3a219bb6, = wb35's binary), 2 x 8 threads. No code change.
# Base = working_branch.env (wb34 stack) + ATM_RH_SIGMA_LAT=52 + ATM_RH_SIGMA_LAT_OCEAN=1 (wb35b; set here, NOT in working_branch.env).
# Control = wb35b: 1033 mm/a (+5.6 %), r .570, sigma 1.36, land 0-15 / 15-35 / 35-65 / 65-90 = 2170 / 9 / 555 / 256
# (NASA 1653 / 643 / 643 / 295), ocean 35-65 1177, 65-90 band 222; S China / SE US / E Austral / SE Africa / S Brazil all 0.00 mm/d
# (NASA 3.87 / 3.89 / 1.86 / 2.22 / 4.55), deserts 0.00, Amazon 8.17, Congo 7.73, max cell 23.9 mm/d.
# STEP 2 BREAKDOWN (wb34 iteration 0): the descent Gaussian (0.74-0.91) is the whole east-side drying; east fetch only 0.67-0.85, so
# EAST 1.0 leaves surface RH 0.73-0.78 (wet end 0.81). At 1.4 the factor (1 - strength*m_east) has NO CLAMP: E Australia -> ~0.86.
# wb26 precedent (older stack, no ceiling): EAST 1.0 gave land 15-35 2 -> 57 and the added moisture rained in the tropics (land 0-15 +320).
# LAND-POLAR criteria (at 600): land 15-35 > 320, land 35-65 515-770, 65-90 band > 180, land 0-15 1400-1900, deserts < 0.3 mm/d, + wb34's.
# PRE-REGISTERED: a (1.0) -- land 15-35 9 -> 30-120, warm east sides (S Brazil, SE Africa, E Austral) 0.1-1.0 mm/d, cool ones (S China, SE US)
# < 0.1, land 0-15 2170 -> 2200-2400 (ceiling limits the spill), deserts < 0.1, global 1035-1055, r .565-.575, ocean unchanged within 1 %;
# b (1.4) -- land 15-35 80-300, warm east sides 0.5-3 mm/d, S China / SE US 0.05-0.8, Arabia up to 0.3, land 0-15 2250-2550, global 1045-1085.
# RISK: Arabia / Oman rain (fetch 0.67); E Australia cells > 10 mm/d at 1.4 (unclamped); land 0-15 moves further above its 1400-1900 range.
set -u; cd "$(dirname "$0")"; rm -f WB36_DONE
for t in wb36a wb36b; do
  mkdir output_$t || { touch WB36_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB36_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_sgl $(md5sum < ../cli/atm_sgl | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_SIGMA_LAT=52 ATM_RH_SIGMA_LAT_OCEAN=1"
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.0 ../cli/atm_sgl config_wb36a.xml > wb36a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.4 ../cli/atm_sgl config_wb36b.xml > wb36b.log 2>&1 &
wait
for t in wb36a wb36b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  echo "banner diff vs wb35b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb35b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean\|P_rain mean\|P_snow mean" $t.log | tail -3
  tail -1 output_$t/convergence.csv
  python3 oceanb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landlat.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
touch WB36_DONE
