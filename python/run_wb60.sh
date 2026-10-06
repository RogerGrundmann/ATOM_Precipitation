#!/bin/bash
# 2026-10-06: wb60a / b / c = working branch (wb57 stack) with a FLATTER storm-track humidity factor (user: "pursue the Southern Ocean 65-90"). Existing knobs:
#   a ATM_RH_STORM 1.15 -> 1.12 + ATM_RH_STORM_POLAR 1.08 -> 1.10;   b 1.10 + 1.10 (flat poleward of 55 deg);   c 1.12 + 1.12.
# SCREENING AT nm 60 from scratch, cli/atm_lev (-O2, md5 02779c82; ATM_LAND_EVAP off = wb57's physics), 3 x 5 threads (wb59 is running on 8).
# WHY (socean.py on wb57, ocean rows S / N against NASA, mm/a): 34-38 deg 424 / 425 (1001 / 1211); 46-50 1392 / 1351 (1033 / 1227); 54-58 1538 / 1980
# (1162 / 1100); 58-62 1279 / 1809 (1086 / 1048); 62-66 552 / 1263 (725 / 988); 66-70 163 / 651 (508 / 783); 70-74 129 / 420 (453 / 444).
# The model's storm-track rain is a steep function of the column water (PW 5-7 mm: 497, 7-10: 1184, 10-14: 2062; NASA 870 / 1114 / 1122): too much at
# 46-62 deg, too little poleward of 62S, and the 35-65 band mean (1184 / 1107) holds by compensation. At the SAME column water both hemispheres rain the
# same in the model, while NASA has twice the rain in the south (PW 2.5-3.5 mm: 487 vs 369; 1.5-2.5: 424 vs 203) -- import by the circumpolar storms,
# which no local quantity knows. So the hemispheric difference is structural; what CAN be tested is the peaked shape.
# Control = wb57: ocean 35-65 1184 (NASA 1107), ocean 65-90 268 (428), Arctic 328 / 379, Southern 184 / 490, global -0.3 %, r .607.
# PRE-REGISTERED (the storm factor is a cliff: 1.0 / 1.15 / 1.25 gave 35-65 25 / 216 / 859 on an older branch): a -- rows 50-62 down 15-30 %, ocean 35-65
# 950-1100, Southern 65-90 200-240, Arctic 380-450, r .605-.612;  b -- rows 50-62 down 30-50 %, ocean 35-65 800-1000, global -2 to -4 %;
# c -- Southern 270-300, Arctic 480-540 (over), ocean 35-65 950-1100.  A usable arm keeps ocean 35-65 within 10 % of NASA and the Arctic below 1.2x NASA.
# RESULT (11:54-12:03, NaN 0; global | r | sigma | ocean 35-65 | land 35-65 | Arctic / Southern ocean 65-90 | S rows 46-50 / 54-58 / 62-66 / 66-70 | N rows 54-58 / 58-62):
#   wb57:           -0.3 % | .607 | 1.23 | 1184 | 664 | 328 / 184 | 1392 / 1538 / 552 / 163 | 1980 / 1809     (NASA 1107 | 643 | 379 / 490 | 1033 / 1162 / 725 / 508 | 1100 / 1048)
#   a 1.12 / 1.10: -10.5 % | .617 | 1.19 |  854 | 483 | 411 / 221 |  971 / 1123 / 486 / 203 | 1396 / 1343
#   b 1.10 / 1.10: -16.1 % | .611 | 1.19 |  683 | 386 | 410 / 220 |  734 /  890 / 484 / 202 | 1076 / 1179
#   c 1.12 / 1.12:  -9.4 % | .612 | 1.20 |  879 | 484 | 520 / 287 |  971 / 1135 / 605 / 269 | 1413 / 1498
# NO ARM IS USABLE (ocean 35-65 must stay within 10 % of NASA): lowering the peak removes the rain the 35-46 deg rows and mid-latitude land already lacked
# (34-38 deg 424 -> 325 against 1000-1200; land 35-65 664 -> 483). But the SHAPE improves: in arm a the southern rows 46-62 deg sit on NASA and r rises
# .607 -> .617, sigma 1.23 -> 1.19. The Southern Ocean 65-90 gains 37 mm/a (a) or 103 with the Arctic 37 % over (c). The factor is centred at 55 deg;
# the observed storm-track rain is flat from 34 to 62 deg -- the centre and width are constants in InitValues_Atm.cpp, not knobs.
set -u; cd "$(dirname "$0")"; rm -f WB60_DONE
for t in wb60a wb60b wb60c; do
  mkdir output_$t || { touch WB60_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB60_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_lev $(md5sum < ../cli/atm_lev | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2"
env OMP_NUM_THREADS=5 $K ATM_RH_STORM=1.12 ATM_RH_STORM_POLAR=1.10 ../cli/atm_lev config_wb60a.xml > wb60a.log 2>&1 &
env OMP_NUM_THREADS=5 $K ATM_RH_STORM=1.10 ATM_RH_STORM_POLAR=1.10 ../cli/atm_lev config_wb60b.xml > wb60b.log 2>&1 &
env OMP_NUM_THREADS=5 $K ATM_RH_STORM=1.12 ATM_RH_STORM_POLAR=1.12 ../cli/atm_lev config_wb60c.xml > wb60c.log 2>&1 &
wait
for t in wb60a wb60b wb60c; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  echo "banner diff vs wb57 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb57.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90
  python3 oceanb.py $V
  python3 socean.py $t 60 | sed -n "2,14p" | cut -c1-200
  python3 polar.py $t 60 | sed -n "16,19p" | cut -c1-110
done
touch WB60_DONE
