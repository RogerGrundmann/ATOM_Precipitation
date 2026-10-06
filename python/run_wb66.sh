#!/bin/bash
# 2026-10-06: wb66 = wb65a AT 600 FROM SCRATCH (user: "a, run the 600, adopt and push when it holds"). Stack = working_branch.env (the wb63 stack, 51 knobs,
# storm factor 1.13 / width 17) + ATM_RH_STORM_SST=0.005, _MAX=0.04, _LAT=48 (_REF at its default 4 C). cli/atm_ss2 (-O2, md5 0fdeda45), 8 threads.
# Control = wb63 (600): 967.6 mm/a (-1.1 %), r .628, sigma 1.18, bands 1769/564/955/262, land/ocean 642/1096, ocean 35-65 1081 (NASA 1107),
# ocean 65-90 264 (428), E Austral 3.20, wettest cell 10.6 mm/d, E 850, P/E 1.14, drift 967.8 -> 967.6.
# PRE-REGISTERED from wb65a at 60: global 962-972, r .628-.633, sigma 1.15-1.17, ocean 35-65 1060-1085, 65-90 band 310-325 (Arctic ocean 395-415,
# Southern 295-315), N rows 54-58 / 58-62 1290-1340 / 1310-1360, bands 0-15 / 15-35, land and every listed region within 0.5 % of wb63, E Austral 3.20,
# drift |100 -> 600| < 1 %, NaN 0.
# RESULT (14:03-14:25, exit 0, NaN 0): 970.1 mm/a (-0.8 %), r .631, sigma 1.16, bands 1769/564/946/316, land/ocean 642.6/1099.7, ocean 35-65 1070, ocean 65-90 364
# (Arctic 405 / 379, Southern 307 / 490), land 65-90 265, N rows 54-58 / 58-62 1301 / 1326, S rows 1228 / 1219, E Austral 3.20, wettest cell 10.6 mm/d,
# E 849.2, P/E 1.14, drift 970.3 -> 970.1. Every pre-registered range met. ADOPTED into working_branch.env.
set -u; cd "$(dirname "$0")"; rm -f WB66_DONE
t=wb66
mkdir output_$t || { touch WB66_DONE; exit 1; }
[ -e config_$t.xml ] && { touch WB66_DONE; exit 1; }
sed -e "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_ss2 $(md5sum < ../cli/atm_ss2 | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_STORM_SST=0.005 ATM_RH_STORM_SST_MAX=0.04 ATM_RH_STORM_SST_LAT=48"
env OMP_NUM_THREADS=8 $K ../cli/atm_ss2 config_$t.xml > $t.log 2>&1
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
echo "banner diff vs wb65a (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb65a.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output" | head -8
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
ITER=520 python3 pacband.py $t | grep -v "^ \+-\?[0-9]\+:\|ocean by latitude"
python3 evapb.py $t 520 | sed -n "2,10p" | cut -c1-130
python3 polar.py $t 520 | sed -n "2,19p" | cut -c1-110
python3 socean.py $t 520 | sed -n "2,14p" | cut -c1-200
python3 polar.py $t 520 | sed -n "2,19p" | cut -c1-110
python3 nrows.py $t 520 | sed -n "/by surface temperature/,/northern ocean 50-66N/p" | cut -c1-110
touch WB66_DONE
