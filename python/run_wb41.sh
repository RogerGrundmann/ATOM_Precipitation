#!/bin/bash
# 2026-10-05: wb41a / wb41b = wb40b stack + ATM_RH_OCEAN_SST=0.007 (/K, ref 29.5 C, largest reduction 0.02), ocean ceiling 0.045 (a) / 0.06 (b).
# User 10-05: "on ocean r-humid is on a high value constant which can't be, consequently the evaporation is low at the equator" and
# "humidity must follow the surface temperature"; NASA's Pacific bands are narrow with centre lines, the model's wide and flat.
# SCREENING: nm 220 from scratch, cli/atm_sst (-O2, 972bbb1 + the knob, md5 3ba679cf), 2 x 8 threads; the -O0 byte check (run_vsst.sh) alongside.
# Stack = working_branch.env + ATM_RH_SIGMA_LAT=52 + _OCEAN=1 + ATM_RH_LAND_EAST=1.2 + ATM_MC_MB_SAT_LAND=0.035 + ATM_HCRIT_SFC=1 + _LAT=10 (wb40b).
# Control = wb40b: 1049 mm/a (+7.3 %), r .573, sigma 1.31, ocean 0-15 2009 (NASA 1440), ocean 15-35 741 (809), land 0-15 1658 (1653), land 15-35 343;
# |lat|<=30 ocean: mean 3.94 (NASA 3.05), r(model, NASA) .655, r(model, SST) .889, p50/p90/p99 4.3/6.8/7.7 (2.5/6.4/8.4), rainless 1 % (7 %);
# by SST 25-26 / 26-27 / 27-28 / 28-29 / 29+ C: 2.31 / 4.20 / 5.74 / 6.49 / 7.00 against NASA 1.77 / 2.32 / 3.71 / 4.79 / 6.55, below 25 C 0.9 vs 1.5-1.7
# (all stratiform); E / C / W Pacific width above half peak 30 / 38 / 43 deg (NASA 9 / 30 / 27), E Pacific 5S-5N 3.9 (1.8); ocean surface RH 0.81 +- 0.01
# over 18S-18N, E 2.7-3.0 mm/d flat. Sensitivity: 0.02 in RH ~ a factor 1.65 in tropical-ocean rain (wb29a vs wb31a).
# PRE-REGISTERED: a -- ocean 0-15 1450-1750, ocean 15-35 480-650 (worse: the cold-water rain is stratiform and RH-sensitive), by SST 26-27 C 2.2-3.0,
# 27-28 C 3.6-4.6, 28-29 C 5.0-5.9, 29+ C unchanged within 3 %; E Pacific width 18-26 deg and 5S-5N 2.0-3.2; r(model, NASA) over |lat|<=30 ocean .66-.70;
# global 900-985, r .575-.60, sigma 1.25-1.40; land 0-15 within 3 % of 1658 (its wet end follows ATM_RH_OCEAN, not the SST); E at 12-18 deg +6-12 %.
# b -- as a with the warm pool allowed to peak: 29+ C 8-10 mm/d, p99 9-11, ocean 0-15 1650-2050, sigma 1.4-1.7, r .56-.60.
# RISK: rain collapses on the cool half of the tropical ocean (rainless area >> NASA's 7 %); the global mean falls below -10 %.
set -u; cd "$(dirname "$0")"; rm -f WB41_DONE
for t in wb41a wb41b; do
  mkdir output_$t || { touch WB41_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB41_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_sst $(md5sum < ../cli/atm_sst | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_SIGMA_LAT=52 ATM_RH_SIGMA_LAT_OCEAN=1"
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.2 ATM_MC_MB_SAT_LAND=0.035 ATM_HCRIT_SFC=1 ATM_HCRIT_SFC_LAT=10 ATM_RH_OCEAN_SST=0.007 ATM_RH_OCEAN_SST_MAX=0.02 ../cli/atm_sst config_wb41a.xml > wb41a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.2 ATM_MC_MB_SAT_LAND=0.035 ATM_HCRIT_SFC=1 ATM_HCRIT_SFC_LAT=10 ATM_RH_OCEAN_SST=0.007 ATM_RH_OCEAN_SST_MAX=0.02 ATM_MC_MB_SAT_OCEAN=0.06 ../cli/atm_sst config_wb41b.xml > wb41b.log 2>&1 &
wait
for t in wb41a wb41b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  echo "banner diff vs wb40b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb40b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean\|P_rain mean\|P_snow mean" $t.log | tail -3
  tail -1 output_$t/convergence.csv
  python3 oceanb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
  python3 landlat.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
python3 pacband.py wb40b wb41a wb41b
touch WB41_DONE
