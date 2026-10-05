#!/bin/bash
# 2026-10-05: wb43a / wb43b = wb40b stack + ATM_RH_OCEAN_SST=0.007, _MAX=0.014, _COLD=25.5 (user: add the cool-side cut-off), ocean ceiling 0.045 (a) / 0.05 (b).
# SCREENING AT nm 60 (user 10-05: the result is there by iteration 20-40; 600 only once at the end): from scratch, cli/atm_ssc (-O2, 2f4e1bb + the
# cut-off, md5 e69847fd), 2 x 8 threads, ATM_VTK_STRIDE=2 for a slice at 60; the -O0 byte check (run_vssc.sh) runs alongside.
# Stack = working_branch.env + ATM_RH_SIGMA_LAT=52 + _OCEAN=1 + ATM_RH_LAND_EAST=1.2 + ATM_MC_MB_SAT_LAND=0.035 + ATM_HCRIT_SFC=1 + _LAT=10 (wb40b).
# wb41a/b (slope 0.007, max 0.02, NO cut-off, ceiling 0.045 / 0.06, at 220): E Pacific width 30 -> 16 deg, 5S-5N 3.9 -> 1.6 / 1.8 (NASA 1.8),
# r(model, NASA) over |lat|<=30 ocean .655 -> .70, BUT rain by SST 25-26 / 26-27 C 2.31 / 4.20 -> 0.16 / 1.71 (NASA 1.77 / 2.32), below 25 C 0.9 -> 0.07,
# rainless 1 % -> 38 % (NASA 7 %), ocean 15-35 741 -> 305 / 369 (809); ceiling 0.06: 29+ C 9.8 vs 6.6, W Pacific peak 10.4 vs 8.4.
# Control = wb40b (220): 1049 mm/a (+7.3 %), r .573, sigma 1.31; by SST <24 / 25-26 / 26-27 / 27-28 / 28-29 / 29+ C: 0.87 / 2.31 / 4.20 / 5.74 / 6.49 / 7.00
# against NASA 1.48 / 1.77 / 2.32 / 3.71 / 4.79 / 6.55.
# PRE-REGISTERED (at 60): a -- by SST <24 C 0.8-0.9, 25-26 C 1.6-2.3, 26-27 C 2.0-2.9, 27-28 C 3.5-4.3, 28-29 C 5.3-6.0; rainless 2-8 %; ocean 15-35 600-720,
# ocean 0-15 1600-1800; E Pacific 5S-5N 2.0-3.0 and width 18-26 deg; r(model, NASA) tropical ocean .67-.70; global 950-1010, r .575-.590, sigma 1.25-1.35;
# b -- 28-29 C 6.2-7.0, 29+ C 7.5-8.3, ocean 0-15 1750-1950, global 985-1040, sigma 1.30-1.45, r .58-.60.
# RISK: the cut-off ramp (25.5 -> 26.5 C) makes a visible rain edge along that isotherm; the E Pacific cold tongue (< 25.5 C) rains again.
set -u; cd "$(dirname "$0")"; rm -f WB43_DONE
for t in wb43a wb43b; do
  mkdir output_$t || { touch WB43_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB43_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_ssc $(md5sum < ../cli/atm_ssc | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2 ATM_RH_SIGMA_LAT=52 ATM_RH_SIGMA_LAT_OCEAN=1 ATM_RH_LAND_EAST=1.2 ATM_MC_MB_SAT_LAND=0.035 ATM_HCRIT_SFC=1 ATM_HCRIT_SFC_LAT=10 ATM_RH_OCEAN_SST=0.007 ATM_RH_OCEAN_SST_MAX=0.014 ATM_RH_OCEAN_SST_COLD=25.5"
env OMP_NUM_THREADS=8 $K ../cli/atm_ssc config_wb43a.xml > wb43a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_MC_MB_SAT_OCEAN=0.05 ../cli/atm_ssc config_wb43b.xml > wb43b.log 2>&1 &
wait
for t in wb43a wb43b; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  echo "banner diff vs wb40b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb40b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -12
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1 | cut -c1-110
  python3 oceanb.py $V
  python3 landb.py $V | grep -v "Oman\|Namib\|W Austral\|Sahara"
done
ITER=60 python3 pacband.py wb40b wb43a wb43b
touch WB43_DONE
