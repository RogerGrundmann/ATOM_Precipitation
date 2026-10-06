#!/bin/bash
# 2026-10-06: wb49 = wb47b AT 600 FROM SCRATCH (user: "adopt wb47b and run the 600, screen a looser cap first"; the screen is run_wb48.sh -- no looser cap
# keeps land 15-35 > 320 with E Australia < 2x NASA). Stack = working_branch.env (wb45) + ATM_RH_LAND_EAST=1.6 + ATM_RH_LAND_EAST_MAX=0.010.
# cli/atm_emx (-O2, md5 3fef9110), 8 threads. Verification run for adopting the two knobs into working_branch.env.
# Control = wb45 (600): 975.4 mm/a (-0.3 %), r .605, sigma 1.26, bands 1754/536/1028/220, land/ocean 667/1098, land 0-15 / 15-35 / 35-65 1627 / 319 / 643,
# E Austral 6.07, wettest cell 18.2 mm/d, no drift (978.2 -> 975.4 over 100-600).
# PRE-REGISTERED from wb47b at 60 (a 60 screen predicted wb45 to 0.3 %): global 962-976, r .605-.615, sigma 1.22-1.26, land 15-35 205-230 (FAILS > 320, known),
# E Austral 2.5-3.0, Madagascar 2.9-3.4, wettest cell < 12 mm/d, no drift (|100 -> 600| < 1 %), NaN 0; every other wb45 criterion kept.
# RESULT (08:16-08:38, exit 0, NaN 0): 968.7 mm/a (-1.0 %), r .609, sigma 1.23, bands 1762/502/1034/220, land/ocean 648/1096, land 0-15 / 15-35 / 35-65
# 1663 / 217 / 664, E Austral 2.70, Madagascar 3.17, S China 0.80, SE US 1.19, S Brazil 0.83, Arabia 0.22, wettest cell 10.5 mm/d, P/E 1.80, converged 1,
# drift 971.4 -> 968.7 over 100-600 (-0.3 %, as wb45). Every pre-registered range met; wb47b at 60 predicted it to 0.3 %. ADOPTED into working_branch.env.
set -u; cd "$(dirname "$0")"; rm -f WB49_DONE
t=wb49
mkdir output_$t || { touch WB49_DONE; exit 1; }
[ -e config_$t.xml ] && { touch WB49_DONE; exit 1; }
sed -e "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_emx $(md5sum < ../cli/atm_emx | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_LAND_EAST=1.6 ATM_RH_LAND_EAST_MAX=0.010"
env OMP_NUM_THREADS=8 $K ../cli/atm_emx config_$t.xml > $t.log 2>&1
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
echo "banner diff vs wb47b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb47b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output" | head -8
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
touch WB49_DONE
