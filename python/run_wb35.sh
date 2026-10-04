#!/bin/bash
# 2026-10-04: wb35a / wb35b = working branch (the wb34 stack) + ATM_RH_SIGMA_LAT=52, land only (a) / land and ocean (b) (user; LAND-POLAR step 1).
# SCREENING: nm 220 from scratch, cli/atm_sgl (-O2, 123e4f0 + the knob), 2 x 8 threads; the -O0 byte check (run_vsgl.sh) runs alongside.
# Control = wb33b / wb34 (same knob set): 971 mm/a (-0.7 %), r .582, sigma 1.36, bands 2036/515/878/37, land/ocean 561/1134;
# land 0-15 / 15-35 / 35-65 / 65-90 = 2165 / 9 / 419 / 8 (NASA 1653 / 643 / 643 / 295), ocean 35-65 1043 (1107), ocean 65-90 62 (428);
# land by latitude (N): 40-52 1.4-2.2 mm/d (NASA 1.5-1.9), 56 0.70, 60 0.25, 64 0.08, 68-76 0.01-0.02 (NASA 1.9 .. 0.6).
# CAUSE: p_sl follows T_s, cold columns stand at ~905 hPa at sea level, MW's sigma = p/p_0 reads them as elevated: surface RH 0.76 / 0.71 at
# 60N / 68N land against H_crit 0.76 / 0.74; with the column's own pressure 0.84 / 0.80.
# LAND-POLAR criteria (user 10-04, at 600): land 15-35 > 320, land 35-65 515-770, 65-90 band > 180, land 0-15 1400-1900, deserts < 0.3 mm/d,
# and the wb34 criteria (global within 10 %, sigma < 2.0, r >= .55, no cell > 50 mm/d).
# PRE-REGISTERED: a -- land 35-65 419 -> 480-600, land 65-90 8 -> 40-150, land 56 / 60 / 64N 0.9-1.6 / 0.8-1.8 / 0.5-1.5 mm/d, land 40-48N unchanged
# within 3 %, ocean and tropics unchanged to the digit, global 978-992, r .58-.59, sigma 1.34-1.36;
# b -- as a on land, plus ocean 35-65 1043 -> 1150-1500, ocean 65-90 62 -> 150-500, global 1010-1120, r .57-.60.
# RISK: elevated land poleward of 52 deg (Altai, northern Rockies, Scandinavia, Greenland, Antarctica) floods -- cells > 10 mm/d;
# b lifts the tuned ocean storm track above NASA.
set -u; cd "$(dirname "$0")"; rm -f WB35_DONE
for t in wb35a wb35b; do
  mkdir output_$t || { touch WB35_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB35_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_sgl $(md5sum < ../cli/atm_sgl | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_SIGMA_LAT=52"
env OMP_NUM_THREADS=8 $K ../cli/atm_sgl config_wb35a.xml > wb35a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_SIGMA_LAT_OCEAN=1 ../cli/atm_sgl config_wb35b.xml > wb35b.log 2>&1 &
wait
for t in wb35a wb35b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  echo "banner diff vs wb34 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb34.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean\|P_rain mean\|P_snow mean" $t.log | tail -3
  tail -1 output_$t/convergence.csv
  python3 oceanb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk | grep -E "Sahara|Amazon|Congo|global max"
  python3 landlat.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
touch WB35_DONE
