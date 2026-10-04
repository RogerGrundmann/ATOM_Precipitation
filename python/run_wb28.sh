#!/bin/bash
# 2026-10-04: wb28a / wb28b = working branch (wb27 stack, now in working_branch.env) + ATM_RAIN_AREA 0.05 / 0.02 (user; STORM screen).
# Default 0.10 (fitted 09-01). Stratiform rain evaporation S_ev scales as f^(5/9): 0.05 -> 0.68x, 0.02 -> 0.41x.
# SCREENING: nm 220 from scratch, cli/atm_g3 (-O2, = cli/atm at 759f496), 2 x 8 threads. Control = wb25b (same knob set, RAIN_AREA 0.10).
# wb25b at 220: 912 mm/a, r .562, sigma 2.85, bands 3001/290/204/4.5, land/ocean 443/1098; ocean 35-65 ~138 (NASA 1107), land 35-65 380 (643).
# Diagnosis (87E, wb27): ocean rain forms at 1.5-4 km and evaporates below ~1.2 km -- 30 % survives at 40-56S, 6 % at 36S; land (1200 m) 96 %.
# PRE-REGISTERED: a (0.05) -- ocean 35-65 138 -> 190-260, 35-65 band 204 -> 245-300; b (0.02) -- ocean 35-65 260-400, band 300-420.
# Both: land 35-65 380 -> 380-420 (no sub-cloud layer on elevated land); 65-90 4.5 -> 5-12; tropics: P_conv unchanged within 1 %
# (convective evaporation is a separate term), 0-15 within 2 %, ocean 15-35 +10-60 mm/a from the 28-36 deg stratiform rain; global a 925-950, b 950-1000;
# r within .01 of .562; sigma falls slightly (2.7-2.85); no NaN. RISK (09-01 record): the response to f was NON-monotonic through accretion and
# snowmelt -- b may not exceed a; and the first use of the new working_branch.env: the banner must equal wb25b's except RAIN_AREA.
set -u; cd "$(dirname "$0")"; rm -f WB28_DONE
for t in wb28a wb28b; do
  mkdir output_$t || { touch WB28_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB28_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_g3 $(md5sum < ../cli/atm_g3 | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1"
env OMP_NUM_THREADS=8 $K ATM_RAIN_AREA=0.05 ../cli/atm_g3 config_wb28a.xml > wb28a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RAIN_AREA=0.02 ../cli/atm_g3 config_wb28b.xml > wb28b.log 2>&1 &
wait
for t in wb28a wb28b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  echo "banner diff vs wb25b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb25b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean\|P_rain mean\|P_snow mean" $t.log | tail -3
  tail -1 output_$t/convergence.csv
  python3 oceanb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk | grep -E "Amazon|Congo|global max"
done
touch WB28_DONE
