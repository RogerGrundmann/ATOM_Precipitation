#!/bin/bash
# 2026-10-05: wb42a / wb42b = wb40b stack with ATM_HCRIT_SFC_LAT 10 -> 25 and ATM_RH_LAND_EAST 1.2 -> 1.25 (a) / 1.35 (b) (user: yes, queue after wb41).
# SCREENING: nm 220 from scratch, cli/atm_sst (-O2, md5 3ba679cf; ATM_RH_OCEAN_SST left at its default 0), 2 x 8 threads. No code change.
# WHY (user 10-05: "Mexico and South Africa have precipitation hotspots"): since wb38 (ATM_RH_LAND_EAST=1.2) the Mexican plateau rains 5.0 mm/d
# (NASA 1.98; wettest cell 19.7 at 24N 100W) and the Highveld 5.5 (1.72; 19.1 at 27S 30E), ALL stratiform (conv 0.00), ground 1540 m: the
# cloud-on-the-ground mechanism outside ATM_HCRIT_SFC_LAT=10 (taper ends at 20 deg). wb40b land 15-35 by elevation, strat | conv | NASA mm/d:
# 0-300 m 0.51 | 0.58 | 2.04;  300-800 0.58 | 0.15 | 1.53;  800-1300 0.73 | 0.03 | 1.50;  1300-2000 2.22 | 0.00 | 1.90 -- the 343 mm/a that
# passed the land 15-35 criterion rests partly on the plateaus (~300 without them); the lowlands, where NASA has most, get least.
# IDEA: limit 25 (full to 25 deg, none from 35) removes plateau cloud at the ground over the subtropics; a stronger east fetch moves the rain
# to the east-coast lowlands (S China 1.15 / 3.87, SE US 0.83 / 3.89).
# Control = wb40b: 1049 mm/a (+7.3 %), r .573, sigma 1.31, land 0-15 1658, land 15-35 343, land 35-65 628, Arabia 0.22, E Austral 3.58, S Brazil 2.62.
# CRITERIA (LAND-POLAR, at 600): land 0-15 1400-1900, land 15-35 > 320, land 35-65 515-770, 65-90 band > 180, deserts < 0.3, global within 10 %, r >= .55.
# PRE-REGISTERED: a (1.25) -- Mexico plt and Highveld < 1.0 mm/d, no land cell > 15 mm/d, land 15-35 230-320, S China 1.4-1.8, SE US 0.9-1.2,
# E Austral 4.2-5.0, S Brazil 2.8-3.8, Arabia 0.2-0.3, land 0-15 1670-1730, land 35-65 560-628 (the taper reaches 35 deg), global 1035-1060, r .570-.585;
# b (1.35) -- land 15-35 360-560, S China 2.4-3.6, SE US 1.5-2.0, E Austral 7-9.5 (over NASA 1.86), S Brazil 6-8.5, Arabia 0.3-0.4, global 1065-1110,
# r .540-.575, wettest land cell 20-35 mm/d.
# RISK: with the plateaus dry, no strength passes land 15-35 > 320 before E Australia / S Brazil flood and r falls below .55 (the 1.4 pattern, wb36b).
set -u; cd "$(dirname "$0")"; rm -f WB42_DONE
for t in wb42a wb42b; do
  mkdir output_$t || { touch WB42_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB42_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_sst $(md5sum < ../cli/atm_sst | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_SIGMA_LAT=52 ATM_RH_SIGMA_LAT_OCEAN=1"
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.25 ATM_MC_MB_SAT_LAND=0.035 ATM_HCRIT_SFC=1 ATM_HCRIT_SFC_LAT=25 ../cli/atm_sst config_wb42a.xml > wb42a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.35 ATM_MC_MB_SAT_LAND=0.035 ATM_HCRIT_SFC=1 ATM_HCRIT_SFC_LAT=25 ../cli/atm_sst config_wb42b.xml > wb42b.log 2>&1 &
wait
for t in wb42a wb42b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  echo "banner diff vs wb40b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb40b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean\|P_rain mean\|P_snow mean" $t.log | tail -3
  tail -1 output_$t/convergence.csv
  python3 oceanb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landlat.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
python3 tropland.py wb40b wb42a wb42b | grep -v elev
touch WB42_DONE
