#!/bin/bash
# 2026-10-05: wb39a / wb39b = wb38b stack + ATM_HCRIT_SFC=1 with ATM_HCRIT_SFC_LAT=15 (a) / every latitude (b) (user: yes, write it; LAND-POLAR step 3).
# SCREENING: nm 220 from scratch, cli/atm_hcs (-O2, 093afd8 + the knob), 2 x 8 threads; the -O0 byte check (run_vhcs.sh) runs alongside.
# Stack = working_branch.env + ATM_RH_SIGMA_LAT=52 + _OCEAN=1 (wb35b) + ATM_RH_LAND_EAST=1.2 + ATM_MC_MB_SAT_LAND=0.035 (wb38b) + the knob.
# Control = wb38b: 1089 mm/a (+11.3 %), r .572, sigma 1.36, land 0-15 2354 (conv 1336, stratiform ~1020; NASA 1653), land 15-35 360,
# land 35-65 628, land 65-90 256, Amazon 7.25, Congo 7.71 (4.75), E Africa 9.96 (conv 2.29; NASA 2.58), Arabia 0.27, max cell 24.4 mm/d.
# CAUSE addressed: hCrit is pressure-only (1.0 at 1000 hPa, 0.64 at 850) and the initial surface RH on a tropical plateau is ~0.69, so cloud
# sits AT elevated ground with no evaporating layer below: tropical land > 800 m (20 % of the area) rains 5-8 mm/d stratiform, no convection.
# CRITERIA (LAND-POLAR, at 600): land 0-15 1400-1900, land 15-35 > 320, land 35-65 515-770, 65-90 band > 180, deserts < 0.3, global within 10 %, r >= .55.
# PRE-REGISTERED: a (lat 15: full to 15 deg, none from 25) -- land 0-15 2354 -> 1750-2050, E Africa 9.96 -> 3-5, Congo 7.71 -> 4.5-6,
# Amazon within 3 %, land 15-35 300-350, global 1060-1075, r .572-.585, everything poleward of 25 deg unchanged to the digit, ocean unchanged;
# b (every latitude) -- as a in the tropics, plus land 15-35 250-330, land 35-65 628 -> 450-580, land 65-90 256 -> 120-230, global 1035-1060.
# RISK: the water no longer rained out near plateau ground rains from a higher deck anyway, or feeds convection on the plateaus (cells > 25 mm/d);
# b undoes step 1 on high-latitude land (Greenland, Antarctica, Siberian uplands are elevated).
set -u; cd "$(dirname "$0")"; rm -f WB39_DONE
for t in wb39a wb39b; do
  mkdir output_$t || { touch WB39_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB39_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_hcs $(md5sum < ../cli/atm_hcs | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_SIGMA_LAT=52 ATM_RH_SIGMA_LAT_OCEAN=1"
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.2 ATM_MC_MB_SAT_LAND=0.035 ATM_HCRIT_SFC=1 ATM_HCRIT_SFC_LAT=15 ../cli/atm_hcs config_wb39a.xml > wb39a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.2 ATM_MC_MB_SAT_LAND=0.035 ATM_HCRIT_SFC=1 ../cli/atm_hcs config_wb39b.xml > wb39b.log 2>&1 &
wait
for t in wb39a wb39b; do
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
python3 tropland.py wb38b wb39a wb39b
touch WB39_DONE
