#!/bin/bash
# *** PROPOSAL, NOT RUN YET (2026-10-08). Dose screen after wb71: the wb71 header below is the evidence. ***
# wb71a / b / c = wb69a stack + ATM_HCRIT_SFC_POLAR=1 + ATM_RH_LAND_ML=1500 with   a strength 0.45;   b 0.5;   c 0.5 and the taper 50..60 deg (both knobs).
# wb71a / b / c = wb69a stack (working branch + ATM_MC_CMB_OCEAN=0.85 + ATM_EVAP_GUST=5) + the two new polar-land knobs (run_vplm.sh byte check):
# SCREENING AT nm 60 from scratch, cli/atm_plm (-O2, HEAD + the knobs), 3 x 6 threads.
# WHY (polarland.py on wb68, land poleward of 60 deg, model / NASA mm/a by ground height): N 0-200 m 106 / 555, 200-500 m 250 / 460, 500-1000 m 743 / 532,
# 1000-2000 m 1174 / 452, > 2000 m 765 / 310; S 91 / 367, 199 / 384, 377 / 375, 525 / 296, 192 / 91. W Siberia 62 / 752, Greenland 661 / 392, E Antarctica 194 / 82.
# Control = wb69a: 981.6 mm/a (+0.3 %), r .646, sigma 1.10, bands 1666/608/1015/320, land 35-65 / 65-90 622 / 269 (NASA 643 / 295).
# EXPECTED (no offline estimate exists; direction only): the threshold lowers ground above ~500 m toward the lowland value, the mixed layer raises ground
# below ~500 m; rain equatorward of 55 deg and all ocean unchanged to 0.5 %. The dose is unknown: the ocean needed 0.4 AND the polar RH floor.
# USABLE if N polar land below 200 m >= 300 and above 1000 m <= 700, Antarctica 150-300, Greenland <= 550, land 35-65 within 560-720, global r >= .645.
# wb70 RESULT (nm 60): both knobs at 0.4 -> N 371 / 381 / 452 / 408 / 219 (NASA 555 / 460 / 532 / 452 / 310), S 253 / 264 / 273 / 302 / 114 (367 / 384 / 375 / 296 / 91),
# r .651, land 35-65 / 65-90 626 / 224; 0.25 -> N 209 / 206 / 256 / 257 / 146. Per 0.15 of strength the lowland gains ~160, ground at 1-2 km ~150 mm/a.
# PRE-REGISTERED by linear extrapolation (N lowland | N 1-2 km | S > 2 km | land 35-65 | land 65-90 | global r):
#   a 0.45: 400-450 | 430-480 | 115-130 | 630-650 | 240-265 | >= .648
#   b 0.5 : 440-510 | 460-540 | 120-140 | 640-670 | 260-300 | >= .645
#   c 0.5, 50..60 deg: as b poleward of 60; land rows 50-60N below 300 m 284 / 388 rise toward 754 / 731; land 35-65 680-760 (NASA 643) -- may overshoot
# USABLE if land 65-90 within 250-330 (NASA 295), E Antarctica <= 160, N Europe <= 1150, land 35-65 within 560-720, global r >= .648, no cell > 12 mm/d.
set -u; cd "$(dirname "$0")"; rm -f WB71_DONE
for t in wb71a wb71b wb71c; do
  mkdir output_$t || { touch WB71_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB71_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_plm $(md5sum < ../cli/atm_plm | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2 ATM_MC_CMB_OCEAN=0.85 ATM_EVAP_GUST=5"
env OMP_NUM_THREADS=6 $K ATM_HCRIT_SFC_POLAR=1 ATM_RH_LAND_ML=1500 ATM_RH_LAND_ML_STRENGTH=0.45 ../cli/atm_plm config_wb71a.xml > wb71a.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_HCRIT_SFC_POLAR=1 ATM_RH_LAND_ML=1500 ATM_RH_LAND_ML_STRENGTH=0.5  ../cli/atm_plm config_wb71b.xml > wb71b.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_HCRIT_SFC_POLAR=1 ATM_HCRIT_SFC_POLAR_LAT=50 ATM_RH_LAND_ML=1500 ATM_RH_LAND_ML_STRENGTH=0.5 ATM_RH_LAND_ML_LAT0=50 ATM_RH_LAND_ML_LAT1=60 ../cli/atm_plm config_wb71c.xml > wb71c.log 2>&1 &
wait
for t in wb71a wb71b wb71c; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  echo "banner diff vs wb69a (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb69a.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "output\|nm=\|VTK_STRIDE" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'water budget closure' $t.log | tail -1 | cut -c1-120
  python3 oceanb.py $V | sed -n '1,4p'
  G=5 python3 polarland.py $t 60 | sed -n '3,44p' | cut -c1-200
  python3 landb.py $V | grep "E Austral\|Amazon\|global max"
done
touch WB71_DONE
