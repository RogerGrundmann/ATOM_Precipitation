#!/bin/bash
# 2026-10-06: wb64a / b / c = the wb62b stack (working branch with ATM_RH_STORM 1.13, width 17 -- the candidate in the wb63 600) + ATM_RH_STORM_SST, a new
# knob (user: "pursue the northern rows 54-62, still 45 % too wet"):   a 0.005 /K, max 0.04;   b 0.007 /K, max 0.05;   c as b with the polar floor 1.08 -> 1.04.
# SCREENING AT nm 60 from scratch, cli/atm_sss (-O2, HEAD 4fb4b40 + the knob), 3 x 5 threads (wb63 is running on 8).
# WHY (nrows.py on wb62b): poleward of 50 deg the model's ocean rain follows the sea temperature and NASA's does not -- model / NASA -6..-2 C 0.26-0.38,
# 1..3 C 0.85-0.98, 3..5 C 1.00-1.31, 5..7 C 1.38, 9..14 C 1.53-1.72, the same in both hemispheres. N Atlantic 54-62N 1868 / 2220 (NASA 1251 / 1296),
# S Pacific 54-58S 1721 (1078), S Atlantic 58-62S 517 (981), S Indian 62-66S 244 (664). The northern rows are simply the warmer rows.
# Control = wb62b: global -1.0 %, r .629, sigma 1.18, ocean 35-65 1090 (1107), ocean 65-90 266 (428; Arctic ~330 / 379, Southern ~185 / 490),
# N rows 54-58 / 58-62 1590 / 1531 (1100 / 1048), S rows 1262 / 1103 (1162 / 1086).
# PRE-REGISTERED (0.01 in RH ~ 11 % of the rain near the peak): a -- N rows 54-62 1300-1450, classes 9-14 C ratio 1.2-1.4 and -6..-2 C 0.5-0.8,
# r .632-.645;  b -- N rows 1200-1350, warm classes 1.05-1.25, cold 0.6-1.0, r .635-.650, Arctic ocean 420-520 (OVER: the cold boost sits on top of the
# polar floor);  c -- Arctic 330-420, Southern ocean 65-90 230-300.  Land and the tropics unchanged; global within 2 % of wb62b.
# USABLE if r >= wb62b's .629, ocean 35-65 within 10 % of NASA, the Arctic below 1.2x NASA, and no ocean row above 1.3x NASA.
# RESULT (13:22-13:31, NaN 0; global | r | ocean 35-65 | Arctic / Southern ocean 65-90 | N rows 46-50 / 54-58 / 58-62 | S rows 46-50 / 54-58 | model / NASA at -6..-2 C and 9..14 C, S and N):
#   wb62b:            -1.0 % | .629 | 1090 | ~330 / ~185 | 1167 / 1590 / 1531 | 1214 / 1262 | 0.38, 0.26 and 1.72, 1.53      (NASA 1107 | 379 / 490 | 1227 / 1100 / 1048 | 1033 / 1162)
#   a 0.005, max .04: -3.7 % | .629 |  965 | 409 / 308 |  766 / 1255 / 1340 |  857 / 1226 | 0.65, 0.43 and 1.13, 0.93
#   b 0.007, max .05: -4.4 % | .625 |  931 | 437 / 355 |  667 / 1149 / 1281 |  743 / 1217 | 0.76, 0.50 and 0.94, 0.75
#   c b + polar 1.04: -4.9 % | .628 |  931 | 305 / 309 |  667 / 1149 / 1281 |  743 / 1217 | as b
# The knob does flatten the rain against the sea temperature, brings the northern rows 54-62 to 1.15-1.28x NASA (from 1.45x), and -- not expected -- lifts
# the Southern Ocean poleward of 65 deg from 185 to 308-355 (NASA 490) with the Arctic at 1.08-1.15x. NOT usable as is: it starts at 40 deg, where the
# rows were already on NASA; they lose 20-40 % (ocean 35-65 -13 to -16 %, global -3.7 to -4.9 %), and r does not move. Next: start it at 48-50 deg (run_wb65.sh).
set -u; cd "$(dirname "$0")"; rm -f WB64_DONE
for t in wb64a wb64b wb64c; do
  mkdir output_$t || { touch WB64_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB64_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_sss $(md5sum < ../cli/atm_sss | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2 ATM_RH_STORM=1.13 ATM_RH_STORM_WIDTH=17"
env OMP_NUM_THREADS=5 $K ATM_RH_STORM_SST=0.005 ATM_RH_STORM_SST_MAX=0.04 ../cli/atm_sss config_wb64a.xml > wb64a.log 2>&1 &
env OMP_NUM_THREADS=5 $K ATM_RH_STORM_SST=0.007 ../cli/atm_sss config_wb64b.xml > wb64b.log 2>&1 &
env OMP_NUM_THREADS=5 $K ATM_RH_STORM_SST=0.007 ATM_RH_STORM_POLAR=1.04 ../cli/atm_sss config_wb64c.xml > wb64c.log 2>&1 &
wait
for t in wb64a wb64b wb64c; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  echo "banner diff vs wb62b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb62b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90
  python3 oceanb.py $V
  python3 socean.py $t 60 | sed -n "2,14p" | cut -c1-200
  python3 nrows.py $t 60 | sed -n "/by surface temperature/,/northern ocean 50-66N/p" | cut -c1-110
  python3 nrows.py $t 60 | grep "5[48]-[56][82]N\|5[48]-[56][82]S" | cut -c1-110
  python3 landb.py $V | grep "E Austral\|S China\|SE US\|Arabia\|Sahara\|Madagascar\|global max"
  python3 polar.py $t 60 | sed -n "16,19p" | cut -c1-110
done
touch WB64_DONE
