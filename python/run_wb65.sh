#!/bin/bash
# 2026-10-06: wb65a / b / c = the wb62b stack + ATM_RH_STORM_SST starting further poleward (new ATM_RH_STORM_SST_LAT; wb64 started at 40 deg and took 20-40 %
# off rows that were on NASA):   a 0.005 /K, max 0.04, start 48;   b the same, start 52;   c 0.007 /K, max 0.05, start 50.   Full weight 10 deg poleward.
# SCREENING AT nm 60 from scratch, cli/atm_ss2 (-O2, HEAD 4fb4b40 + the knobs), 3 x 5 threads.
# Control = wb62b: global -1.0 %, r .629, sigma 1.18, ocean 35-65 1090 (NASA 1107), N rows 46-50 / 50-54 / 54-58 / 58-62 1167 / 1387 / 1590 / 1531
# (1227 / 1202 / 1100 / 1048), S rows 1214 / 1243 / 1262 / 1103 (1033 / 1052 / 1162 / 1086), Arctic ~330 / 379, Southern ocean 65-90 ~185 / 490.
# wb64a (start 40): N rows 766 / 967 / 1255 / 1340, Arctic 409, Southern 308, ocean 35-65 965, r .629.
# PRE-REGISTERED: a -- rows 46-50 within 5 % of wb62b, rows 54-62 N 1250-1400, ocean 35-65 1020-1090, global -1.5 to -2.5 %, Arctic 400-415, Southern
# ocean 65-90 300-315, r .632-.645;  b -- rows up to 54 untouched, rows 58-62 N 1350-1450, ocean 35-65 1050-1100;  c -- N rows 54-62 1150-1300, Arctic 430-445.
# USABLE if r >= .629, ocean 35-65 within 10 % of NASA, the Arctic below 1.2x NASA (455), no ocean row above 1.3x NASA.
# RESULT (13:32-13:41, NaN 0; global | r | sigma | ocean 35-65 | 65-90 band | Arctic / Southern ocean 65-90 | N rows 46-50 / 50-54 / 54-58 / 58-62 | S rows 54-58 / 62-66):
#   wb62b:                -1.0 % | .629 | 1.18 | 1090 | 266 | ~330 / ~185 | 1167 / 1387 / 1590 / 1531 | 1262 / 505     (NASA 1107 | 364 | 379 / 490 | 1227 / 1202 / 1100 / 1048 | 1162 / 725)
#   a .005, .04, start 48: -0.8 % | .631 | 1.16 | 1079 | 320 | 409 / 308 | 1159 / 1258 / 1316 / 1341 | 1242 / 689
#   b .005, .04, start 52: -0.3 % | .629 | 1.17 | 1097 | 320 | 409 / 308 | 1167 / 1380 / 1486 / 1383 | 1262 / 688
#   c .007, .05, start 50: -0.3 % | .630 | 1.16 | 1092 | 339 | 437 / 355 | 1166 / 1329 / 1328 / 1297 | 1252 / 764
# a and c are USABLE on every criterion (b leaves the 54-58N row at 1.35x NASA). Rows equatorward of 50 deg, land, the tropics and E Australia (3.20)
# are untouched. The gain is in the high latitudes: 65-90 band 266 -> 320 / 339 (NASA 364), Southern ocean 185 -> 308 / 355, northern rows 54-62 from
# 1.45x to 1.20-1.28x NASA; r +.001-.002, sigma 1.18 -> 1.16.
set -u; cd "$(dirname "$0")"; rm -f WB65_DONE
for t in wb65a wb65b wb65c; do
  mkdir output_$t || { touch WB65_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB65_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_ss2 $(md5sum < ../cli/atm_ss2 | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2 ATM_RH_STORM=1.13 ATM_RH_STORM_WIDTH=17"
env OMP_NUM_THREADS=5 $K ATM_RH_STORM_SST=0.005 ATM_RH_STORM_SST_MAX=0.04 ATM_RH_STORM_SST_LAT=48 ../cli/atm_ss2 config_wb65a.xml > wb65a.log 2>&1 &
env OMP_NUM_THREADS=5 $K ATM_RH_STORM_SST=0.005 ATM_RH_STORM_SST_MAX=0.04 ATM_RH_STORM_SST_LAT=52 ../cli/atm_ss2 config_wb65b.xml > wb65b.log 2>&1 &
env OMP_NUM_THREADS=5 $K ATM_RH_STORM_SST=0.007 ATM_RH_STORM_SST_LAT=50 ../cli/atm_ss2 config_wb65c.xml > wb65c.log 2>&1 &
wait
for t in wb65a wb65b wb65c; do
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
touch WB65_DONE
