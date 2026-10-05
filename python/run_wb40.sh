#!/bin/bash
# 2026-10-05: wb40a / wb40b = wb38b stack + ATM_HCRIT_SFC=1 (SECOND form: ground >= 1000 hPa left alone) with ATM_HCRIT_SFC_LAT=15 (a) / 10 (b).
# SCREENING: nm 220 from scratch, cli/atm_hc2 (-O2, 093afd8 + the knob, md5 031c4394), 2 x 8 threads; the -O0 byte check (run_vhc2.sh) runs alongside.
# wb39a/b (FIRST form, lat 15 / every latitude): highland fixed -- tropical land > 800 m stratiform 5.9-8.0 -> 0.00-0.03 mm/d, E Africa 9.96 -> 3.16
# (NASA 2.58), Congo 7.71 -> 4.79 (4.75), global max leaves the E African highlands -- BUT lowland stratiform ROSE, < 150 m 1.38 -> 5.29 mm/d
# (Amazon 7.25 -> 9.81, Maritime 6.69 -> 9.86), because a warm lowland column stands at ~1040 hPa and 1000/p_g lowers its threshold aloft;
# land 0-15 only 2354 -> 2082. Without the latitude limit (b) extratropical land rain vanishes (land 35-65 628 -> 18, land 65-90 256 -> 2):
# the step-1 high-latitude land rain IS cloud at the ground of cold low-pressure columns. The limit is mandatory.
# Control = wb38b: 1089 mm/a (+11.3 %), r .572, sigma 1.36, land 0-15 2354, land 15-35 360, Amazon 7.25, Congo 7.71, E Africa 9.96, max cell 24.4.
# CRITERIA (LAND-POLAR, at 600): land 0-15 1400-1900, land 15-35 > 320, land 35-65 515-770, 65-90 band > 180, deserts < 0.3, global within 10 %, r >= .55.
# PRE-REGISTERED: a (lat 15) -- tropical land stratiform < 400 m as wb38b (1.4-1.5 mm/d), > 800 m as wb39a (< 0.1); Amazon 7.2-7.4, Maritime 6.6-6.8,
# Congo 4.8-5.5, E Africa 3.1-3.4, land 0-15 1850-1980, land 15-35 310-320, global 1058-1066, r .578-.590, poleward of 25 deg unchanged to the digit;
# b (lat 10: full to 10 deg, none from 20) -- land 0-15 1900-2050, land 15-35 340-360, SE Africa 2.4-2.5, otherwise as a.
set -u; cd "$(dirname "$0")"; rm -f WB40_DONE
for t in wb40a wb40b; do
  mkdir output_$t || { touch WB40_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB40_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_hc2 $(md5sum < ../cli/atm_hc2 | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_SIGMA_LAT=52 ATM_RH_SIGMA_LAT_OCEAN=1"
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.2 ATM_MC_MB_SAT_LAND=0.035 ATM_HCRIT_SFC=1 ATM_HCRIT_SFC_LAT=15 ../cli/atm_hc2 config_wb40a.xml > wb40a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.2 ATM_MC_MB_SAT_LAND=0.035 ATM_HCRIT_SFC=1 ATM_HCRIT_SFC_LAT=10 ../cli/atm_hc2 config_wb40b.xml > wb40b.log 2>&1 &
wait
for t in wb40a wb40b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  echo "banner diff vs wb38b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb38b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean\|P_rain mean\|P_snow mean" $t.log | tail -3
  tail -1 output_$t/convergence.csv
  python3 oceanb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landlat.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
python3 tropland.py wb38b wb39a wb40a wb40b
touch WB40_DONE
