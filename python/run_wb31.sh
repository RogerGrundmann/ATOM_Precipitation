#!/bin/bash
# 2026-10-04: wb31a / wb31b = working branch (wb27 stack) + ATM_RH_OCEAN=0.82 + the SAME mass-flux ceiling on ocean and land,
# ATM_MC_MB_SAT_OCEAN = ATM_MC_MB_SAT_LAND = 0.045 / 0.035 (user; RAIN-CONV). SCREENING: nm 220 from scratch, cli/atm_mbs (-O2, 2f1af19), 2 x 8 threads.
# Control = wb30a (RH_OCEAN 0.82, ocean ceiling 0.065, no land ceiling): 1030 mm/a, r .537, sigma 2.27, bands 3020/617/234/4.5, land/ocean 816/1115,
# ocean 0-15 2883 (conv 2483), ocean 15-35 860, land 0-15 3507, Amazon 16.2, Congo 8.7 mm/d; |lat|<=30 ocean mean 5.39, p50/p90/p99 5.2/10.7/12.0
# (NASA 3.05; 2.6/6.5/8.4); mean M_u(base) 0.033; rain ~185 mm/d per unit M_b at this humidity.
# CLOSING CRITERIA (at 600): global within 10 % of NASA, ocean 0-15 < 2200, sigma < 2.0, r >= .55, no cell > 50 mm/d.
# PRE-REGISTERED: the mean uncapped M_b (~0.037) is BELOW both ceilings, so the mean falls little and only the top is cut:
# a (0.045) -- ocean<=30 mean 4.6-5.1, p99 8.0-8.8, p50 4.5-5.2, ocean 0-15 2550-2800, global 880-960, sigma 1.5-1.8;
# b (0.035) -- ocean<=30 mean 4.1-4.6, p99 6.2-6.9, p50 4.0-4.8, ocean 0-15 2150-2450, global 800-880, sigma 1.3-1.6.
# Land now capped: Amazon 16.2 -> a 7.5-9, b 6-7 mm/d; Congo 8.7 -> a 6.5-8, b 5.5-6.5; land 0-15 3507 -> a 2000-2500, b 1700-2100; no cell > 12 mm/d.
# Both: ocean 15-35 700-860, zero-rain area 0 %, r .53-.56, deserts 0, 35-65 / 65-90 as wb30a, no NaN.
# EXPECTED SHORTFALL: the median stays ~1.6-2x NASA -- the humidity then has to come down (0.81) at the chosen ceiling.
# RISK: the land gate (|M_u| > 0.01) with mean land M_b 0.0085 -- wb22 (blend 1) lost all land convection under a land ceiling.
set -u; cd "$(dirname "$0")"; rm -f WB31_DONE
for t in wb31a wb31b; do
  mkdir output_$t || { touch WB31_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB31_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_mbs $(md5sum < ../cli/atm_mbs | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.82"
env OMP_NUM_THREADS=8 $K ATM_MC_MB_SAT_OCEAN=0.045 ATM_MC_MB_SAT_LAND=0.045 ../cli/atm_mbs config_wb31a.xml > wb31a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_MC_MB_SAT_OCEAN=0.035 ATM_MC_MB_SAT_LAND=0.035 ../cli/atm_mbs config_wb31b.xml > wb31b.log 2>&1 &
wait
for t in wb31a wb31b; do
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
touch WB31_DONE
