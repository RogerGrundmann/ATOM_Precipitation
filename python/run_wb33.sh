#!/bin/bash
# 2026-10-04: wb33a / wb33b = the wb31a stack + ATM_RH_OCEAN_ML=1500 at ATM_RH_OCEAN_ML_STRENGTH 0.25 / 0.4 (user; STORM).
# SCREENING: nm 220 from scratch, cli/atm_oms (-O2, d8193cd + the strength knob), 2 x 8 threads; the -O0 byte check (run_vmoms.sh) runs alongside.
# wb31a (layer off): ocean 35-65 171 (NASA 1107), ocean 65-90 2 (428), land 35-65 401 (643), global 750, r .514, sigma 1.40.
# wb32b (1500 m, fully mixed): ocean 35-65 9732, ocean 65-90 1340, global 3179, r .233; 40-56S: RH - H_crit +0.13 / +0.22 / +0.14 at 600 / 1000 /
# 1500 m, q_c 0.13 / 0.35 / 0.06 g/kg, 31 mm/d at the sea. Baseline RH there 0.70 / 0.67 / 0.61 against H_crit 0.75 / 0.68 / 0.56.
# STORM pass criteria (proposed): ocean 35-65 > 500, land 35-65 <= 770, r >= .55.
# PRE-REGISTERED (q_c ~ (RH - H_crit)^2): a (0.25) -- RH - H_crit -0.005 / +0.045 / +0.07: cloud down to ~900 m, 40-56S surface rain 1.2-2.5 mm/d,
# ocean 35-65 450-900, ocean 65-90 20-120, global 820-920, r .53-.57; b (0.4) -- RH - H_crit +0.02 / +0.08 / +0.09: cloud down to ~550 m,
# 40-56S 3-5 mm/d, ocean 35-65 1000-2200, ocean 65-90 100-400, global 950-1250, r .50-.56. Both: land 35-65 405-440, tropics unchanged to the digit,
# ocean 15-35 709 -> 710-740, no NaN. RISK: accretion makes the response steeper than quadratic (b overshoots again).
set -u; cd "$(dirname "$0")"; rm -f WB33_DONE
for t in wb33a wb33b; do
  mkdir output_$t || { touch WB33_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB33_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_oms $(md5sum < ../cli/atm_oms | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.82 ATM_MC_MB_SAT_OCEAN=0.045 ATM_MC_MB_SAT_LAND=0.045 ATM_RH_OCEAN_ML=1500"
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN_ML_STRENGTH=0.25 ../cli/atm_oms config_wb33a.xml > wb33a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN_ML_STRENGTH=0.4 ../cli/atm_oms config_wb33b.xml > wb33b.log 2>&1 &
wait
for t in wb33a wb33b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  echo "banner diff vs wb31a (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb31a.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean\|P_rain mean" $t.log | tail -2
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-RG\] ocean<15" $t.log | tail -1 | cut -c1-260
  tail -1 output_$t/convergence.csv
  python3 oceanb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk | grep -E "India|Sahara|Arabia|Amazon|Congo|global max"
  python3 tropdist.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk | head -1
  python3 storm87.py output_$t/0Ma_smooth_Atm_zonal_87_220.vtk 2>/dev/null | grep -A8 "OCEAN 40-56S" | cut -c1-150
done
touch WB33_DONE
