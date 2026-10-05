#!/bin/bash
# 2026-10-05: wb38a / wb38b = wb35b base + ATM_RH_LAND_EAST=1.2 (user: carry 1.2, go to step 3), land mass-flux ceiling 0.045 (a) / 0.035 (b).
# LAND-POLAR step 3 (tropical land), first screen: nm 220 from scratch, cli/atm_sgl (-O2, md5 3a219bb6), 2 x 8 threads. No code change.
# a measures 1.2 itself (wb37: 1.15 / 1.25 -> land 15-35 292 / 467, Arabia 0.25 / 0.30, r .578 / .568, global +12.0 / +14.1 %, land 0-15 2602 / 2653).
# READ-ONLY BREAKDOWN (tropland.py + elevation, wb35b / wb37a at 220): tropical land (<15 deg) by ground elevation, strat | conv | NASA mm/d:
#   wb37a  0-150 m 1.49 | 7.13 | 6.29;  150-400 1.42 | 5.21 | 4.67;  400-800 2.73 | 3.87 | 3.73;  800-1300 5.76 | 0.83 | 3.44;  1300-2000 7.84 | 0.05 | 3.24
# -> TWO causes: lowland (54 % of the area) is too CONVECTIVE (+0.5..0.8 mm/d over NASA before the stratiform part); highland (> 800 m, 20 %)
# does not convect and rains 5-8 mm/d STRATIFORM straight onto the ground (E Africa 10.3 vs 2.58 mm/d, Congo strat 4.3). The land ceiling can
# only reach the first; wb31 (0.045 -> 0.035, ocean and land together): land 0-15 2170 -> 1941, Amazon 8.17 -> 6.88, Congo 7.73 -> 7.42.
# CRITERIA (LAND-POLAR, at 600): land 0-15 1400-1900, land 15-35 > 320, deserts < 0.3 mm/d, global within 10 %, r >= .55.
# PRE-REGISTERED: a -- land 15-35 350-400, land 0-15 2620-2640, Arabia 0.26-0.29, global 1103-1110, r .572-.578;
# b -- land 0-15 2350-2450 (still > 1900: the highland stratiform part stays), Amazon 7.2-7.6, Congo 7.6-7.9, E Africa unchanged within 3 %,
# land 15-35 330-390, global 1085-1098, r .570-.580, ocean bands unchanged within 1 %.
set -u; cd "$(dirname "$0")"; rm -f WB38_DONE
for t in wb38a wb38b; do
  mkdir output_$t || { touch WB38_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB38_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_sgl $(md5sum < ../cli/atm_sgl | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_SIGMA_LAT=52 ATM_RH_SIGMA_LAT_OCEAN=1"
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.2 ../cli/atm_sgl config_wb38a.xml > wb38a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.2 ATM_MC_MB_SAT_LAND=0.035 ../cli/atm_sgl config_wb38b.xml > wb38b.log 2>&1 &
wait
for t in wb38a wb38b; do
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
python3 tropland.py wb35b wb38a wb38b
touch WB38_DONE
