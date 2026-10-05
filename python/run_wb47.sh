#!/bin/bash
# 2026-10-05: wb47a / wb47b = working branch (the wb45 stack) + ATM_RH_LAND_EAST_MAX (user: "option 2" = clamp at a land cap slightly above the wet end,
# then retune): a = cap +0.015 at strength 1.35, b = cap +0.010 at strength 1.6.
# SCREENING AT nm 60 from scratch, cli/atm_emx (-O2, 12ed221 + the knob), 2 x 8 threads, ATM_VTK_STRIDE=2; the -O0 byte check (run_vemx.sh) alongside.
# WHY (eastfetch.py on wb45): the factor (1 - 1.35*m_east) is negative where the fetch is nearly complete, so the initial land RH (E Australia 0.846,
# max 0.908) exceeds the wet end 0.82: E Australia 6.07 vs NASA 1.86 mm/d (5.04 stratiform), Madagascar 7.3 vs 4.3; 16 % of the 15-35 deg land area
# carries 87 % of the band's rain (4.80 mm/d there vs NASA 3.39, the rest 0.13 vs 1.46).
# Control = wb45 (600; = wb44a at 60): 975.4 mm/a (-0.3 %), r .605, sigma 1.26, land 0-15 / 15-35 / 35-65 1627 / 319 / 643, E Austral 6.07, S China 1.84,
# SE US 1.63, S Brazil 1.69, Madagascar 7.30, Arabia 0.17, wettest cell 18.2 mm/d (23S 149E).
# ESTIMATE from the run's own rain-vs-RH response (lowland < 500 m at 15-35 deg, mean now 1.42 mm/d; the curve overstates Arabia ~3x and S China):
#   a 1.35 +0.015: lowland 0.99, E Austral 3.2, Madagascar 4.0;   b 1.6 +0.010: lowland 1.20, E Austral 3.2, S China higher, Arabia up ~1.8x.
# PRE-REGISTERED (at 60): a -- E Austral 2.8-3.8, Madagascar 3.5-4.8, no land cell > 12 mm/d at 15-35 deg, land 15-35 319 -> 200-250, S China / SE US /
# S Brazil within 10 % of now, Arabia 0.15-0.19, global 955-968, r .605-.615, sigma 1.22-1.26, land 0-15 within 2 %;
# b -- land 15-35 250-320, E Austral 2.8-3.8, S China 2.2-3.0, SE US 1.8-2.3, S Brazil 2.0-2.8, Arabia 0.25-0.40 (may break the desert limit 0.3),
# global 962-980, r .60-.615.
# RISK: the cliff -- rain at RH 0.83-0.835 may be far from the interpolated 3-4 mm/d; b puts rain on Arabia.
# RESULT (14:20-14:25, NaN 0): a 964 mm/a (-1.4 %), r .609, sigma 1.24, land 15-35 185, E Austral 2.86, Madagascar 3.56, S China 0.65, SE US 0.96, S Brazil 0.50,
# Mexico plt 0.21, Arabia 0.16, wettest cell 10.7 mm/d; b 972 (-0.7 %), r .609, sigma 1.24, land 15-35 217, E Austral 2.73, Madagascar 3.16, S China 0.81,
# SE US 1.19, S Brazil 0.85, Arabia 0.21. The overshoot is removed (wettest cell 18.4 -> 10.7) but land 15-35 falls well below 320: S China, S Brazil and
# the Mexican plateau rain was overshoot cells too. The response-curve estimate was too high for every east-side region.
set -u; cd "$(dirname "$0")"; rm -f WB47_DONE
for t in wb47a wb47b; do
  mkdir output_$t || { touch WB47_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB47_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_emx $(md5sum < ../cli/atm_emx | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2"
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST_MAX=0.015 ../cli/atm_emx config_wb47a.xml > wb47a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST_MAX=0.010 ATM_RH_LAND_EAST=1.6 ../cli/atm_emx config_wb47b.xml > wb47b.log 2>&1 &
wait
for t in wb47a wb47b; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  echo "banner diff vs wb45 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb45.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90
  python3 oceanb.py $V
  python3 landb.py $V | grep -v "Oman\|Namib\|land \|output"
done
touch WB47_DONE
