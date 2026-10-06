#!/bin/bash
# 2026-10-06: wb59 = wb58 AT 600 FROM SCRATCH (user: "run the 600, adopt and push when it holds"). Stack = working_branch.env (the wb57 stack, which
# already carries ATM_EVAP_WIND / _GUST and ATM_RH_STORM_POLAR) + ATM_LAND_EVAP=1. cli/atm_lev (-O2, md5 02779c82), 8 threads.
# Control = wb57 (600): 975.0 mm/a (-0.3 %), r .607, sigma 1.23, bands 1766/504/1035/261, land/ocean 648/1104, E 758.6 (land 0), P/E 1.29,
# E Austral 2.70, wettest cell 10.5 mm/d, drift 975.8 -> 975.0 over 100-600.
# PRE-REGISTERED: land E 315-335 mm/a (E/P 0.48-0.52), global E 845-860 (ocean E falls ~2 % over the run, as in wb54 / wb57), P/E 1.13-1.16;
# rain within 0.3 % of wb57 in the global mean, every band, land / ocean and every listed region; r .607; same drift; NaN 0.
# RESULT (11:51-12:21, slowed by the wb60 screen; exit 0, NaN 0): 975.5 mm/a (-0.3 %), r .607, sigma 1.23, bands 1766.5/504.6/1035.9/261.4, land/ocean 649.7/1104.5,
# land E 325.2 (E/P 0.50, potential 1125), ocean E 1062, global E 850.8, P/E 1.15, E Austral 2.71, wettest cell 10.6 mm/d, converged 1, drift 975.9 -> 975.5.
# Every pre-registered range met. ADOPTED into working_branch.env.
set -u; cd "$(dirname "$0")"; rm -f WB59_DONE
t=wb59
mkdir output_$t || { touch WB59_DONE; exit 1; }
[ -e config_$t.xml ] && { touch WB59_DONE; exit 1; }
sed -e "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_lev $(md5sum < ../cli/atm_lev | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_LAND_EVAP=1"
env OMP_NUM_THREADS=8 $K ../cli/atm_lev config_$t.xml > $t.log 2>&1
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
echo "banner diff vs wb58 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb58.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output" | head -8
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
touch WB59_DONE
