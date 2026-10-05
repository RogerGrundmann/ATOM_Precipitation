#!/bin/bash
# 2026-10-05: wb44a / wb44b = the best land stack (wb42b) + the best ocean stack (wb43b), ATM_RH_LAND_EAST 1.35 (a) / 1.3 (b). Candidate for the 600.
# SCREENING AT nm 60 from scratch, cli/atm_ssc (-O2, md5 e69847fd), 2 x 8 threads, ATM_VTK_STRIDE=2. No code change.
# Stack = working_branch.env + ATM_RH_SIGMA_LAT=52 + _OCEAN=1 + ATM_MC_MB_SAT_LAND=0.035 + ATM_HCRIT_SFC=1 + ATM_HCRIT_SFC_LAT=25 + ATM_RH_LAND_EAST
# + ATM_RH_OCEAN_SST=0.007 + _MAX=0.014 + _COLD=25.5 + ATM_MC_MB_SAT_OCEAN=0.05.
# wb42b (land, 220): +7.7 %, r .577, sigma 1.31, land 0-15 1691, land 15-35 339, land 35-65 645, Arabia 0.25, Mexico plt 1.14, Highveld 1.23,
# E Austral 6.22 (NASA 1.86; wettest cell 18.4), S China 1.85, SE US 1.64. wb43b (ocean, 60): -0.4 %, r .600, sigma 1.26, ocean 1099,
# ocean 0-15 1787, ocean 15-35 614, land 0-15 1597, Arabia 0.15 (the SST knob costs land 0-15 ~4 % and Arabia ~0.07).
# CRITERIA (LAND-POLAR + wb34, at 600): land 0-15 1400-1900, land 15-35 > 320, land 35-65 515-770, 65-90 band > 180, deserts < 0.3 mm/d,
# global within 10 %, r >= .55, sigma < 2.0, ocean 0-15 < 2200, ocean 35-65 > 500, no cell > 50 mm/d.
# PRE-REGISTERED (at 60): a -- global 965-995, r .595-.610, sigma 1.24-1.30, land 0-15 1600-1650, land 15-35 300-340, land 35-65 635-650, Arabia 0.16-0.22,
# E Austral 5.5-6.3, Mexico plt 1.0-1.2, ocean as wb43b within 1 %; b -- land 15-35 260-300 (below the criterion), E Austral 4.6-5.3, S China 1.3-1.6, r as a.
set -u; cd "$(dirname "$0")"; rm -f WB44_DONE
for t in wb44a wb44b; do
  mkdir output_$t || { touch WB44_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB44_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_ssc $(md5sum < ../cli/atm_ssc | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2 ATM_RH_SIGMA_LAT=52 ATM_RH_SIGMA_LAT_OCEAN=1 ATM_MC_MB_SAT_LAND=0.035 ATM_HCRIT_SFC=1 ATM_HCRIT_SFC_LAT=25 ATM_RH_OCEAN_SST=0.007 ATM_RH_OCEAN_SST_MAX=0.014 ATM_RH_OCEAN_SST_COLD=25.5 ATM_MC_MB_SAT_OCEAN=0.05"
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.35 ../cli/atm_ssc config_wb44a.xml > wb44a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST=1.3 ../cli/atm_ssc config_wb44b.xml > wb44b.log 2>&1 &
wait
for t in wb44a wb44b; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  echo "banner diff vs wb40b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb40b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -12
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1 | cut -c1-110
  python3 oceanb.py $V
  python3 landb.py $V | grep -v "Oman\|Namib\|W Austral\|Sahara"
done
ITER=60 python3 pacband.py wb44a wb44b | grep -v "^ \+-\?[0-9]\+:\|ocean by latitude"
touch WB44_DONE
