#!/bin/bash
# 2026-10-07: wb68 = wb67b AT 600 FROM SCRATCH (user: "run b at 600 from scratch"). Stack = working_branch.env (the wb66 stack, 54 knobs) +
# ATM_RH_OCEAN_ML_LAT0=22, ATM_RH_OCEAN_ML_LAT1=38 (shipped 30 / 40). cli/atm_mlt (-O2, md5 44f6333f), 8 threads.
# Control = wb66 (600): 970.1 mm/a (-0.8 %), r .631, sigma 1.16, bands 1769/564/946/316, land/ocean 644/1099, ocean 15-35 693 (NASA 809), ocean 35-65
# 1070 (1107), S ocean rows 22-26 / 26-30 / 30-34 / 34-38 / 38-42 424 / 443 / 462 / 564 / 1074, E Austral 3.20, wettest cell 10.6 mm/d, P/E 1.14.
# PRE-REGISTERED from wb67b at 60: global 1012-1023 (+3.5 to +4.5 %), r .645-.651, sigma 1.13-1.15, ocean 15-35 795-815, ocean 35-65 1150-1172,
# S rows 26-30 / 30-34 / 34-38 / 38-42 540-565 / 815-845 / 1050-1090 / 1175-1200, N row 22-26 705-730, bands 0-15 and 65-90, land and every listed
# land region within 0.5 % of wb66 except E Austral 3.25-3.31, drift |100 -> 600| < 1 %, NaN 0.
# RESULT (09:33-09:55, exit 0, NaN 0): 1016.9 mm/a (+4.0 %), r .648, sigma 1.14, bands 1769/649/1006/316, land/ocean 643.9/1164.6, ocean 15-35 810
# (NASA 809), ocean 35-65 1152 (1107), S rows 22-26 / 26-30 / 30-34 / 34-38 / 38-42 442 / 549 / 827 / 1064 / 1180 (NASA 728 / 841 / 950 / 1001 / 1077),
# N rows 728 / 663 / 826 / 1052 / 1183 (637 / 796 / 1009 / 1211 / 1378), E Austral 3.26, wettest cell 10.6 mm/d, E 848.6, P/E 1.20 (wb66 1.14),
# drift 1017.3 -> 1016.9, converged 1. Every pre-registered range met. At 87E the rain reaching the sea at 28 / 32 / 36S: 23 / 43 / 61 % (wb66 19 / 22 / 43).
# NOT adopted by this script: the user decides after the result.
set -u; cd "$(dirname "$0")"; rm -f WB68_DONE
t=wb68
mkdir output_$t || { touch WB68_DONE; exit 1; }
[ -e config_$t.xml ] && { touch WB68_DONE; exit 1; }
sed -e "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_mlt $(md5sum < ../cli/atm_mlt | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN_ML_LAT0=22 ATM_RH_OCEAN_ML_LAT1=38"
env OMP_NUM_THREADS=8 $K ../cli/atm_mlt config_$t.xml > $t.log 2>&1
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
echo "banner diff vs wb66 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb66.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output" | head -8
echo "-- trajectory (every 100: bands; score; land / ocean)"
grep -a "by |latitude|" $t.log | grep -v MFC | awk 'NR%100==0' | cut -c1-150
grep -a "model .*NASA .*bias" $t.log | awk 'NR%100==0' | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | awk 'NR%100==0' | cut -c1-90
echo "-- final"
grep -a "by |latitude|" $t.log | grep -v MFC | tail -2; grep -a "model .*NASA .*bias" $t.log | tail -2 | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1 | cut -c1-110; grep -a "LAND-EVAP" $t.log | tail -1
grep -a "P_conv mean\|P_rain mean" $t.log | tail -2
grep -a "max v-component\|max w-component\|max u-component" $t.log | tail -3 | cut -c1-170
tail -1 output_$t/convergence.csv
for it in 120 520; do V=output_$t/0Ma_smooth_Atm_radial_0_$it.vtk; echo "-- slice $it"; python3 oceanb.py $V; python3 landb.py $V | grep -v "land \|output"; done
python3 row36.py $t 520 | grep -v "T [0-9]" | cut -c1-230
ITER=520 python3 pacband.py $t | grep -v "^ \+-\?[0-9]\+:\|ocean by latitude"
python3 evapb.py $t 520 | sed -n "2,10p" | cut -c1-130
python3 polar.py $t 520 | sed -n "2,19p" | cut -c1-110
touch WB68_DONE
