#!/bin/bash
# 2026-10-05: wb37a / wb37b = wb35b base + ATM_RH_LAND_EAST 1.15 / 1.25 (user approved 10-05; LAND-POLAR step 2, fills the wb36 gap).
# SCREENING: nm 220 from scratch, cli/atm_sgl (-O2, md5 3a219bb6), 2 x 8 threads. No code change. Base as run_wb36.sh.
# wb36 (1.0 / 1.4): land 15-35 138 / 877 (NASA 643; control wb35b 9), of which convective 71 / 153; land 0-15 2529 / 2736 (2170);
# global 1075 / 1162 (+9.9 / +18.7 %), r .578 / .513, sigma 1.38 / 1.62; S China 0.17 / 3.91, SE US 0.23 / 2.00, E Austral 0.93 / 9.90,
# SE Africa 0.74 / 6.79, S Brazil 0.43 / 8.99, Arabia 0.14 / 0.43 mm/d; max cell 24.5 / 46.7 mm/d. A cliff between 1.0 and 1.4.
# CRITERIA watched here (LAND-POLAR, at 600): land 15-35 > 320, r >= .55, deserts (Arabia) < 0.3 mm/d, global within 10 %, no cell > 50 mm/d.
# PRE-REGISTERED (convex response assumed): a (1.15) -- land 15-35 220-420, global 1090-1115, r .555-.575, Arabia 0.18-0.28, E Austral 2-5 mm/d,
# S China 0.5-1.5, land 0-15 2580-2650; b (1.25) -- land 15-35 400-650, global 1110-1140 (> +10 %), r .535-.560, Arabia 0.25-0.35,
# E Austral 4-7, S Brazil 3-6, S China 1.5-2.8, land 0-15 2630-2700, max cell 28-40 mm/d.
# RISK: no value passes land 15-35 > 320 with global within 10 % (the global mean is already +9.9 % at 1.0) -> needs the step-3 / RC-RESID trim.
set -u; cd "$(dirname "$0")"; rm -f WB37_DONE
for t in wb37a wb37b; do
  mkdir output_$t || { touch WB37_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB37_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_sgl $(md5sum < ../cli/atm_sgl | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_SIGMA_LAT=52 ATM_RH_SIGMA_LAT_OCEAN=1"
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.15 ../cli/atm_sgl config_wb37a.xml > wb37a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.25 ../cli/atm_sgl config_wb37b.xml > wb37b.log 2>&1 &
wait
for t in wb37a wb37b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  echo "banner diff vs wb35b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb35b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean\|P_rain mean\|P_snow mean" $t.log | tail -3
  tail -1 output_$t/convergence.csv
  python3 oceanb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landlat.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
touch WB37_DONE
