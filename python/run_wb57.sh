#!/bin/bash
# 2026-10-06: wb57 = ONE 600 FROM SCRATCH on the stack to keep (user: "1.08, combine with the evaporation knobs into one 600"):
# working_branch.env (wb49) + ATM_EVAP_WIND=1 + ATM_EVAP_GUST=3 (wb53b) + ATM_RH_STORM_POLAR=1.08 (wb56b). cli/atm_spl (-O2, md5 71233269), 8 threads.
# Control = wb49 (600): 968.7 mm/a (-1.0 %), r .609, sigma 1.23, bands 1762/502/1034/220, land/ocean 648/1096, land 0-15 / 15-35 / 35-65 1663 / 217 / 664,
# E 537.5, P/E 1.80, E Austral 2.70, wettest cell 10.5 mm/d, drift 971.4 -> 968.7 over 100-600.
# PRE-REGISTERED (the two screens are independent: E is a null on the rain, the polar factor acts on ocean poleward of ~67 deg only): global 970-976,
# r .606-.609, sigma 1.23-1.24, 65-90 band 255-270 (ocean 260-275, Arctic 320-340, Southern 178-190), bands 0-15 / 15-35 / 35-65 within 0.5 % of wb49,
# land and every listed region within 1 %, E 770-790, P/E 1.23-1.27, wettest cell < 12 mm/d, drift |100 -> 600| < 1 %, NaN 0.
# RESULT (11:04-11:30, exit 0, NaN 0): 975.0 mm/a (-0.3 %), r .607, sigma 1.23, bands 1766/504/1035/261, land/ocean 648/1104, land 0-15 / 15-35 / 35-65
# 1663 / 217 / 664, ocean 65-90 268 (Arctic 328 / 379, Southern 184 / 490), land 65-90 255, E 758.6 (ocean 1062, land 0), P/E 1.29, E Austral 2.70,
# wettest cell 10.5 mm/d, converged 1, drift 975.8 -> 975.0 over 100-600. Every pre-registered range met except E / P-E (770-790 / 1.23-1.27; E falls
# slightly over the run, as in wb54). ADOPTED into working_branch.env.
set -u; cd "$(dirname "$0")"; rm -f WB57_DONE
t=wb57
mkdir output_$t || { touch WB57_DONE; exit 1; }
[ -e config_$t.xml ] && { touch WB57_DONE; exit 1; }
sed -e "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_spl $(md5sum < ../cli/atm_spl | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_EVAP_WIND=1 ATM_EVAP_GUST=3 ATM_RH_STORM_POLAR=1.08"
env OMP_NUM_THREADS=8 $K ../cli/atm_spl config_$t.xml > $t.log 2>&1
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
echo "banner diff vs wb56b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb56b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output" | head -8
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
python3 evapb.py $t 520 | sed -n "2,10p" | cut -c1-130
python3 polar.py $t 520 | sed -n "2,19p" | cut -c1-110
touch WB57_DONE
