#!/bin/bash
# 2026-10-06: wb63 = wb62b AT 600 FROM SCRATCH (user: "b, run the 600, adopt and push when it holds"). Stack = working_branch.env (the wb59 stack, 50 knobs,
# incl. ATM_LAND_EVAP=1) with ATM_RH_STORM 1.15 -> 1.13 and ATM_RH_STORM_WIDTH 15 -> 17. cli/atm_sto (-O2, md5 b76e4aa1), 8 threads.
# Control = wb59 (600): 975.5 mm/a (-0.3 %), r .607, sigma 1.23, bands 1767/505/1036/261, land/ocean 650/1105, land 0-15 / 15-35 / 35-65 1663 / 217 / 664,
# ocean 35-65 1184 (NASA 1107), ocean 15-35 614 (809), E 851, P/E 1.15, E Austral 2.71, wettest cell 10.6 mm/d, drift 975.9 -> 975.5.
# PRE-REGISTERED from wb62b at 60: global 958-972 (-1 to -2 %), r .625-.632, sigma 1.17-1.19, ocean 35-65 1070-1100, land 35-65 610-630, ocean 15-35
# 680-695, land 15-35 235-245, 65-90 band 258-268, E Austral 3.1-3.3 (KNOWN COST), SE US 1.5-1.6, wettest cell < 12 mm/d, E 845-860, P/E 1.12-1.15,
# drift |100 -> 600| < 1 %, NaN 0.
# RESULT (13:20-14:01, slowed by the wb64 / wb65 screens; exit 0, NaN 0): 967.6 mm/a (-1.1 %), r .628 (.629 at 100-500), sigma 1.18, bands 1769/564/955/262,
# land/ocean 642.3/1096.4, ocean 15-35 / 35-65 693 / 1081, land 15-35 / 35-65 240 / 619, E Austral 3.20, SE US 1.51, S China 0.93, wettest cell 10.6 mm/d,
# E 850.2, P/E 1.14, converged 1, drift 967.8 -> 967.6 over 100-600. Every pre-registered range met. ADOPTED into working_branch.env.
set -u; cd "$(dirname "$0")"; rm -f WB63_DONE
t=wb63
mkdir output_$t || { touch WB63_DONE; exit 1; }
[ -e config_$t.xml ] && { touch WB63_DONE; exit 1; }
sed -e "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_sto $(md5sum < ../cli/atm_sto | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_STORM=1.13 ATM_RH_STORM_WIDTH=17"
env OMP_NUM_THREADS=8 $K ../cli/atm_sto config_$t.xml > $t.log 2>&1
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
echo "banner diff vs wb62b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb62b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output" | head -8
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
touch WB63_DONE
