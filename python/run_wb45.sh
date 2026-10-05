#!/bin/bash
# 2026-10-05: wb45 = wb44a to 600 from scratch (user: "run wb44a to 600") -- the verification of the best stack against the LAND-POLAR and wb34 criteria.
# cli/atm_ssc (-O2, 0fe4e96, md5 e69847fd), 8 threads, default VTK stride (slices at 20, 120, ... 520).
# Stack = working_branch.env (wb34) + ATM_RH_SIGMA_LAT=52 + _OCEAN=1 + ATM_RH_LAND_EAST=1.35 + ATM_MC_MB_SAT_LAND=0.035 + ATM_HCRIT_SFC=1 + ATM_HCRIT_SFC_LAT=25
# + ATM_RH_OCEAN_SST=0.007 + _MAX=0.014 + _COLD=25.5 + ATM_MC_MB_SAT_OCEAN=0.05.
# wb44a at 60: 978.5 mm/a (+0.0 %), r .605, sigma 1.26, bands 1752/534/1039/223, land/ocean 670/1101, land 0-15 / 15-35 / 35-65 / 65-90 1628 / 321 / 647 / 257,
# ocean 0-15 1787, ocean 35-65 1185, Arabia 0.16, E Austral 6.15, wettest cell 18.4 mm/d, P/E 1.79.
# CRITERIA (at 600): global within 10 %, r >= .55, sigma < 2.0, ocean 0-15 < 2200, ocean 35-65 > 500, no cell > 50 mm/d (wb34);
# land 15-35 > 320, land 35-65 515-770, 65-90 band > 180, land 0-15 1400-1900, Sahara / Arabia / Namib / W Australia < 0.3 mm/d (LAND-POLAR).
# PRE-REGISTERED: every number within 1 % of its value at 60 (wb34 moved 972.3 -> 970.0 over 100-600; wb40b 1052 -> 1049 over 40-220): global 968-980,
# r .600-.608, sigma 1.25-1.27, land 15-35 316-324 -- MARGINAL against > 320 --, max |v| ~2.6 m/s (no seam mode), converged 1, NaN 0.
# RISK: land 15-35 ends below 320; a slow drift appears in 0-15 as on wb27 (+0.7 % over 600).
# RESULT (12:10-12:32, 22 min, exit 0, NaN 0, converged 1, max |v| 2.6 m/s): 975.4 mm/a (-0.3 %), r .605, sigma 1.26, bands 1754/536/1028/220, land/ocean 667/1098,
# P/E 1.81; no drift (global 978.2 -> 975.4 over 100-600, all of it 35-65 1038 -> 1028 and 65-90 222.9 -> 219.8, slowing); at 520: land 0-15 1627, land 15-35 319,
# land 35-65 643, ocean 0-15 1788, ocean 35-65 1173, Arabia 0.17, E Austral 6.07, wettest cell 18.2 mm/d. EVERY criterion met EXCEPT land 15-35 (319 vs > 320).
set -u; cd "$(dirname "$0")"; rm -f WB45_DONE
t=wb45
mkdir output_$t || { touch WB45_DONE; exit 1; }
[ -e config_$t.xml ] && { touch WB45_DONE; exit 1; }
sed -e "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_ssc $(md5sum < ../cli/atm_ssc | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_SIGMA_LAT=52 ATM_RH_SIGMA_LAT_OCEAN=1 ATM_RH_LAND_EAST=1.35 ATM_MC_MB_SAT_LAND=0.035 ATM_HCRIT_SFC=1 ATM_HCRIT_SFC_LAT=25 ATM_RH_OCEAN_SST=0.007 ATM_RH_OCEAN_SST_MAX=0.014 ATM_RH_OCEAN_SST_COLD=25.5 ATM_MC_MB_SAT_OCEAN=0.05"
env OMP_NUM_THREADS=8 $K ../cli/atm_ssc config_$t.xml > $t.log 2>&1
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
echo "banner diff vs wb44a (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb44a.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output" | head -8
echo "-- trajectory (every 100: bands; score; land / ocean)"
grep -a "by |latitude|" $t.log | grep -v MFC | awk 'NR%100==0' | cut -c1-150
grep -a "model .*NASA .*bias" $t.log | awk 'NR%100==0' | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | awk 'NR%100==0' | cut -c1-90
echo "-- final"
grep -a "by |latitude|" $t.log | grep -v MFC | tail -2; grep -a "model .*NASA .*bias" $t.log | tail -2 | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1 | cut -c1-110
grep -a "P_conv mean\|P_rain mean" $t.log | tail -2
grep -a "max v-component\|max w-component\|max u-component" $t.log | tail -3 | cut -c1-170
tail -1 output_$t/convergence.csv
for it in 120 520; do V=output_$t/0Ma_smooth_Atm_radial_0_$it.vtk; echo "-- slice $it"; python3 oceanb.py $V; python3 landb.py $V | grep -v "land \|output"; done
ITER=520 python3 pacband.py $t | grep -v "^ \+-\?[0-9]\+:\|ocean by latitude"
touch WB45_DONE
