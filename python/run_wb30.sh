#!/bin/bash
# 2026-10-04: wb30a / wb30b = working branch (wb27 stack) + ATM_MC_MB_SAT_OCEAN=0.065 with ATM_RH_OCEAN 0.82 / 0.85 (user; RAIN-CONV).
# SCREENING: nm 220 from scratch, cli/atm_mbs (-O2, 2f1af19), 2 x 8 threads. Controls: wb25b (no ceiling, RH_OCEAN 0.80), wb29a (0.045, 0.80).
# wb25b: 912 mm/a, r .562, sigma 2.85, ocean 0-15 / 15-35 3385 / 406, ocean 1098, land 443, land 0-15 1641, Amazon 8.25, Congo 3.12.
# wb29a: 454 mm/a, r .524, sigma 1.12, ocean 0-15 / 15-35 1212 / 255, ocean 458; ceiling measured at 120 mm/d per unit M_s -> 0.065 ~ 7.8 mm/d.
# |lat|<=30 ocean p50/p75/p90/p99: NASA 2.6/4.6/6.5/8.4, wb29a 1.9/4.1/5.0/5.4; 25 % rains 0.
# CLOSING CRITERIA for RAIN-CONV (user 10-04, at 600): global within 10 % of NASA, ocean 0-15 < 2200, sigma < 2.0, r >= .55, no cell > 50 mm/d.
# PRE-REGISTERED: a (0.82) -- ocean<=30 mean 3.6-4.4 mm/d, p99 7.5-8.0, zero-rain area 25 % -> 14-20 %, ocean 0-15 1900-2300, ocean 15-35 450-650,
# ocean 720-880, global 700-850, sigma 1.6-2.0, r .54-.58; b (0.85) -- ocean<=30 mean 4.4-5.3, zero-rain area 5-12 %, ocean 0-15 2200-2600,
# ocean 15-35 650-950, ocean 880-1080, global 850-1050, sigma 1.8-2.2, r .54-.58. LAND IS NOT CAPPED and its wet end follows ATM_RH_OCEAN
# (ATM_RH_LAND=2): land 0-15 1641 -> a 2000-2600, b 3000-4500; Amazon a 10-13, b 14-20 mm/d; land cells > 50 mm/d possible in b. Deserts 0;
# 35-65, 65-90 unchanged (RH_OCEAN acts within 30 deg). No NaN. RISK: land overshoot dominates sigma; the humidity response is still steep.
set -u; cd "$(dirname "$0")"; rm -f WB30_DONE
for t in wb30a wb30b; do
  mkdir output_$t || { touch WB30_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB30_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_mbs $(md5sum < ../cli/atm_mbs | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_MC_MB_SAT_OCEAN=0.065"
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN=0.82 ../cli/atm_mbs config_wb30a.xml > wb30a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN=0.85 ../cli/atm_mbs config_wb30b.xml > wb30b.log 2>&1 &
wait
for t in wb30a wb30b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  echo "banner diff vs wb25b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb25b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean\|P_rain mean" $t.log | tail -2
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-RG\] ocean<15" $t.log | tail -1 | cut -c1-260
  tail -1 output_$t/convergence.csv
  python3 oceanb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk | grep -E "India|Sahara|Arabia|Amazon|Congo|global max"
  python3 tropdist.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
touch WB30_DONE
