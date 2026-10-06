#!/bin/bash
# 2026-10-06: wb51a / b / c = working branch (wb49 stack) + ATM_RH_LAND_EAST_ML=1500, refining wb50 (0.3 too little, 0.6 too much):
#   a strength 0.4;  b strength 0.5;  c strength 0.5 with ATM_RH_LAND_EAST_ML_T=23 (full weight at <= 21 C instead of <= 22 C).
# SCREENING AT nm 60 from scratch, cli/atm_eml (-O2, md5 b933515b), 3 x 6 threads. Control = wb49 / wb50 (see run_wb50.sh).
# PRE-REGISTERED (interpolating wb50a/b; the response steepens): a -- S China 2.3-3.2, SE US 2.3-3.0, S Brazil 1.8-2.5, E Austral 5.0-6.0, land 15-35
# 340-400, r .603-.609;  b -- S China 3.2-4.3, SE US 3.0-3.6, S Brazil 2.4-3.0, E Austral 6.0-7.2, land 15-35 400-460, r .598-.606;
# c -- E Austral (22.5 C) and S Brazil (22.3 C) well below b (weight x0.2-0.4), S China / SE US as b, Highveld below b.
# RESULT (09:48-09:54, NaN 0; land 15-35 / S China / SE US / S Brazil / E Austral / Highveld / Mexico plt / global / r / sigma / wettest cell):
#   a 0.4:        350 / 2.48 / 2.44 / 1.93 / 5.26 / 1.94 / 0.80 / +0.8 % / .608 / 1.25 / 12.6
#   b 0.5:        416 / 3.48 / 3.11 / 2.55 / 6.45 / 2.63 / 1.17 / +1.5 % / .603 / 1.26 / 15.8
#   c 0.5, T 23:  388 / 3.18 / 3.06 / 2.06 / 5.57 / 2.58 / 1.07 / +1.2 % / .604 / 1.26 / 15.8
# All inside the pre-registered ranges. Land 15-35 passes 320 from 0.4 with r .608. The cooler threshold separates little: E Australia -14 %, S Brazil -19 %.
# E Australia is the price at every strength (2.3x / 2.8x / 3.5x NASA at 0.3 / 0.4 / 0.5).
set -u; cd "$(dirname "$0")"; rm -f WB51_DONE
for t in wb51a wb51b wb51c; do
  mkdir output_$t || { touch WB51_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB51_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_eml $(md5sum < ../cli/atm_eml | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2 ATM_RH_LAND_EAST_ML=1500"
env OMP_NUM_THREADS=6 $K ATM_RH_LAND_EAST_ML_STRENGTH=0.4 ../cli/atm_eml config_wb51a.xml > wb51a.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_RH_LAND_EAST_ML_STRENGTH=0.5 ../cli/atm_eml config_wb51b.xml > wb51b.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_RH_LAND_EAST_ML_STRENGTH=0.5 ATM_RH_LAND_EAST_ML_T=23 ../cli/atm_eml config_wb51c.xml > wb51c.log 2>&1 &
wait
for t in wb51a wb51b wb51c; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  grep -a "RH-LAND-EAST-ML" $t.log | head -1
  echo "banner diff vs wb49 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb49.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90
  python3 oceanb.py $V
  python3 landb.py $V | grep -v "Oman\|Namib\|land \|output"
done
touch WB51_DONE
