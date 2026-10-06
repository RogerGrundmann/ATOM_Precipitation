#!/bin/bash
# 2026-10-06: wb52 = wb51a AT 600 FROM SCRATCH (user: "0.4, run the 600, adopt and push when it holds"). Stack = working_branch.env (wb49) +
# ATM_RH_LAND_EAST_ML=1500 + ATM_RH_LAND_EAST_ML_STRENGTH=0.4 (T at its default 24). cli/atm_eml (-O2, md5 b933515b), 8 threads.
# Control = wb49 (600): 968.7 mm/a (-1.0 %), r .609, sigma 1.23, bands 1762/502/1034/220, land/ocean 648/1096, land 0-15 / 15-35 / 35-65 1663 / 217 / 664,
# E Austral 2.70, wettest cell 10.5 mm/d, drift 971.4 -> 968.7 over 100-600.
# PRE-REGISTERED from wb51a at 60: global 978-990, r .603-.612, sigma 1.23-1.27, land 15-35 335-365, land 35-65 675-700, S China 2.3-2.6, SE US 2.3-2.6,
# S Brazil 1.8-2.1, E Austral 5.0-5.5 (KNOWN COST, 2.8x NASA), Highveld 1.8-2.1, wettest cell < 14 mm/d, no drift (|100 -> 600| < 1 %), NaN 0.
# STOPPED BY THE USER at iteration ~160 (10:11, killed by PID; output_wb52 kept, incomplete): "wb49 was much better". At 160: 985.3 mm/a (+0.7 %), r .608,
# sigma 1.25 = wb51a. ATM_RH_LAND_EAST_ML is NOT adopted; the working branch stays the wb49 stack.
set -u; cd "$(dirname "$0")"; rm -f WB52_DONE
t=wb52
mkdir output_$t || { touch WB52_DONE; exit 1; }
[ -e config_$t.xml ] && { touch WB52_DONE; exit 1; }
sed -e "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_eml $(md5sum < ../cli/atm_eml | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_LAND_EAST_ML=1500 ATM_RH_LAND_EAST_ML_STRENGTH=0.4"
env OMP_NUM_THREADS=8 $K ../cli/atm_eml config_$t.xml > $t.log 2>&1
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
echo "banner diff vs wb51a (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb51a.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output" | head -8
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
touch WB52_DONE
