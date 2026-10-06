#!/bin/bash
# 2026-10-06: wb54 = wb53b AT 600 FROM SCRATCH (user: "b, run the 600, adopt and push when it holds"). Stack = working_branch.env (wb49) +
# ATM_EVAP_WIND=1 + ATM_EVAP_GUST=3. cli/atm_evw (-O2, md5 b88a0004), 8 threads.
# Control = wb49 (600): 968.7 mm/a (-1.0 %), r .609, sigma 1.23, bands 1762/502/1034/220, land/ocean 648/1096, land 0-15 / 15-35 / 35-65 1663 / 217 / 664,
# E 537.5, P/E 1.80, E Austral 2.70, wettest cell 10.5 mm/d, drift 971.4 -> 968.7 over 100-600.
# PRE-REGISTERED: E 765-785, P/E 1.23-1.27; rain within 0.5 % of wb49 in the global mean, every band, land / ocean and every listed region, r .609,
# sigma 1.23 (E adds ~0.0008 mm more in 120 s); same drift as wb49; NaN 0.
# RESULT (10:35-11:11, slowed by the wb55/wb56 screens; exit 0, NaN 0): E 759.3 mm/a (ocean 1063, land 0), P/E 1.28 -- just outside the pre-registered
# 765-785 / 1.23-1.27 (E falls slightly as level 1 moistens: 776 at 60). Rain 971.0 mm/a (-0.7 %), r .609, sigma 1.24, E Austral 2.70, wettest cell 10.5:
# as wb49 to 0.3 %, and the drift over 100-600 is smaller (971.7 -> 971.0; wb49 971.4 -> 968.7). Superseded for adoption by wb57 (the combined stack).
set -u; cd "$(dirname "$0")"; rm -f WB54_DONE
t=wb54
mkdir output_$t || { touch WB54_DONE; exit 1; }
[ -e config_$t.xml ] && { touch WB54_DONE; exit 1; }
sed -e "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_evw $(md5sum < ../cli/atm_evw | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_EVAP_WIND=1 ATM_EVAP_GUST=3"
env OMP_NUM_THREADS=8 $K ../cli/atm_evw config_$t.xml > $t.log 2>&1
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
echo "banner diff vs wb53b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb53b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output" | head -8
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
ITER=520 python3 evapb.py $t 520 | sed -n "2,10p" | cut -c1-130
touch WB54_DONE
