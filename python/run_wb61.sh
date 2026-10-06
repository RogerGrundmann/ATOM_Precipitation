#!/bin/bash
# 2026-10-06: wb61a / b / c = working branch (wb57 stack) with the storm-track humidity factor MOVED EQUATORWARD (new knobs; user: "pursue the storm-track
# shape"). The factor stays 1.15 and the polar floor 1.08:   a ATM_RH_STORM_LAT 55 -> 48;   b -> 45;   c -> 45 with ATM_RH_STORM_WIDTH 15 -> 12.
# SCREENING AT nm 60 from scratch, cli/atm_sto (-O2, HEAD 3d93007 + ATM_LAND_EVAP (off) + the two knobs), 3 x 6 threads, after wb59 has finished.
# WHY (socean.py on wb57; wb60): ocean rain S / N against NASA, mm/a -- 34-38 deg 424 / 425 (1001 / 1211), 38-42 951 / 944 (1077 / 1378), 46-50 1392 / 1351
# (1033 / 1227), 54-58 1538 / 1980 (1162 / 1100), 58-62 1279 / 1809 (1086 / 1048), 62-66 552 / 1263 (725 / 988). The factor peaks at 55 deg; the observed
# rain is flat from 34 to 62. Lowering the peak (wb60) fixed the shape (r .617) and lost the mean. Factor at 36 / 45 / 55 / 62 deg: now 1.030 / 1.096 /
# 1.150 / 1.121; centre 48: 1.079 / 1.144 / 1.121 / 1.080 (floor); centre 45: 1.105 / 1.150 / 1.096 / 1.080; centre 45, width 12: 1.085 / 1.150 / 1.080 / 1.080.
# It also reaches the subtropics: at 27 deg 1.005 now, 1.021 (48), 1.036 (45), 1.016 (45 / 12) -- E Australia and ocean 15-35 will feel it.
# Control = wb57: global -0.3 %, r .607, sigma 1.23, ocean 35-65 1184 (NASA 1107), land 35-65 664 (643), ocean 15-35 614 (809), E Austral 2.70 (1.86).
# PRE-REGISTERED (a cliff, so wide): rows 34-42 deg up by 30-100 %, rows 54-62 down by 15-40 %, ocean 35-65 within 10 % of NASA in a and c;
# r .610-.625; E Austral a 3.0-3.6, b 3.3-4.5, c 2.9-3.4; ocean 15-35 up 3-12 %; land 35-65 550-700.
# USABLE if ocean 35-65 and land 35-65 stay within 10 % of NASA, r >= .612 and E Australia <= 3.2 mm/d.
# RESULT (12:21-12:28, NaN 0; global | r | sigma | ocean 35-65 | land 35-65 | ocean 15-35 | S rows 34-38 / 38-42 / 46-50 / 54-58 | E Austral):
#   wb57:              -0.3 % | .607 | 1.23 | 1184 |  664 |  614 |  424 /  951 / 1392 / 1538 | 2.70      (NASA 1107 | 643 | 809 | 1001 / 1077 / 1033 / 1162 | 1.86)
#   a centre 48:      +34.6 % | .539 | 1.47 | 1983 | 1115 | 1032 | 1887 / 3038 / 2203 / 1105 | 5.33
#   b centre 45:      +60.0 % | .440 | 1.85 | 2375 | 1326 | 1547 | 3226 / 4282 / 2121 /  828 | 8.10
#   c centre 45, w 12: +34.1 % | .489 | 1.65 | 2058 | 1153 |  932 | 2142 / 3685 / 2053 /  712 | 4.82
# ALL THREE FAIL, far outside the pre-registered ranges: a factor of 1.08-1.15 on the warm, water-rich columns at 34-46 deg rains 2000-4500 mm/a.
# The response per row gives the factor each latitude needs (now -> needed): 36 deg 1.030 -> ~1.045, 40 deg 1.055 -> ~1.06, 48 deg 1.121 -> ~1.11,
# 56 deg 1.149 -> ~1.13, 60 deg 1.134 -> ~1.12: not a shift but a slightly LOWER and WIDER factor (run_wb62.sh).
set -u; cd "$(dirname "$0")"; rm -f WB61_DONE
for t in wb61a wb61b wb61c; do
  mkdir output_$t || { touch WB61_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB61_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_sto $(md5sum < ../cli/atm_sto | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2"
env OMP_NUM_THREADS=6 $K ATM_RH_STORM_LAT=48 ../cli/atm_sto config_wb61a.xml > wb61a.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_RH_STORM_LAT=45 ../cli/atm_sto config_wb61b.xml > wb61b.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_RH_STORM_LAT=45 ATM_RH_STORM_WIDTH=12 ../cli/atm_sto config_wb61c.xml > wb61c.log 2>&1 &
wait
for t in wb61a wb61b wb61c; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  echo "banner diff vs wb57 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb57.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90
  python3 oceanb.py $V
  python3 socean.py $t 60 | sed -n "2,14p" | cut -c1-200
  python3 landb.py $V | grep "E Austral\|S China\|SE US\|Arabia\|Sahara\|Madagascar\|global max"
  python3 polar.py $t 60 | sed -n "16,19p" | cut -c1-110
done
touch WB61_DONE
