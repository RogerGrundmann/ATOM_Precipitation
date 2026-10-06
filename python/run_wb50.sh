#!/bin/bash
# 2026-10-06: wb50a / b / c = working branch (the wb49 stack) + ATM_RH_LAND_EAST_ML=1500 at strength 0.3 / 0.6 / 1.0 (new knob; user: "option A, write
# the knob and screen it"). SCREENING AT nm 60 from scratch, cli/atm_eml (-O2, HEAD ad5c68c + the knob), 3 x 6 threads, ATM_VTK_STRIDE=2; -O0 byte check
# (run_veml.sh) alongside.
# WHY (eastside.py on wb49): east-side land at the cap rains by its ground temperature -- >= 24 C convects and is on NASA, < 24 C has zero convective
# rain and 1.2-1.3 mm/d stratiform against 3.0-3.3. The knob gives that cool land a (partly) mixed boundary layer, weight strength x m_east.
# Control = wb49 (600; = wb47b at 60): 968.7 mm/a (-1.0 %), r .609, sigma 1.23, land 0-15 / 15-35 / 35-65 1663 / 217 / 664, E Austral 2.70 (NASA 1.86),
# S China 0.80 (3.87), SE US 1.19 (3.89), S Brazil 0.83 (4.55), SE Africa 1.16 (2.22), Madagascar 3.17, Mexico plt 0.21 (1.98), Arabia 0.22, wettest 10.5.
# PRE-REGISTERED (at 60; wide -- the ocean analogue went from too little at 0.25 to 3-9x too much fully mixed):
#   a 0.3: S China / SE US / S Brazil x1.3-2, E Austral 3.0-4.0, land 15-35 240-330;  b 0.6: those three x2-4, E Austral 4-7, land 15-35 300-550;
#   c 1.0: overshoot (some region > 3x NASA, a cell > 20 mm/d), r below .60.  Unchanged everywhere: land 0-15 within 2 %, land 35-65 within 5 %,
#   Arabia / Sahara / W Australia (hot or no fetch), ocean bands.
# QUESTION: gain of S China + SE US + S Brazil per unit of E Australia overshoot, and is there a strength with land 15-35 > 320 and r >= .60?
# RESULT (09:41-09:47, NaN 0; land 15-35 / S China / SE US / S Brazil / E Austral / Highveld / Mexico plt / global / r / sigma / wettest cell):
#   wb49 (off): 217 / 0.80 / 1.19 / 0.83 / 2.70 / 0.51 / 0.21 / -1.0 % / .609 / 1.23 / 10.5
#   a 0.3:      299 / 1.77 / 1.97 / 1.43 / 4.33 / 1.41 / 0.55 / +0.3 % / .609 / 1.24 / 10.8
#   b 0.6:      501 / 4.87 / 3.94 / 3.35 / 7.88 / 3.52 / 1.68 / +2.4 % / .594 / 1.29 / 19.7
#   c 1.0:     1052 / 16.0 / 9.07 / 8.80 / 14.2 / 9.66 / 5.38 / +7.9 % / .482 / 1.64 / 48.6
# The knob does what it was written for: 0.3 doubles S China, r unchanged; 0.6 puts SE US on NASA and S China 26 % above. Land 0-15, the ocean bands,
# Arabia untouched; land 35-65 664 -> 682 / 702 / 740. E Australia rises with the rest (as built), Highveld doubles NASA at 0.6. Next: 0.4-0.5 and a
# cooler threshold (run_wb51.sh).
set -u; cd "$(dirname "$0")"; rm -f WB50_DONE
for t in wb50a wb50b wb50c; do
  mkdir output_$t || { touch WB50_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB50_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_eml $(md5sum < ../cli/atm_eml | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2 ATM_RH_LAND_EAST_ML=1500"
env OMP_NUM_THREADS=6 $K ATM_RH_LAND_EAST_ML_STRENGTH=0.3 ../cli/atm_eml config_wb50a.xml > wb50a.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_RH_LAND_EAST_ML_STRENGTH=0.6 ../cli/atm_eml config_wb50b.xml > wb50b.log 2>&1 &
env OMP_NUM_THREADS=6 $K ATM_RH_LAND_EAST_ML_STRENGTH=1.0 ../cli/atm_eml config_wb50c.xml > wb50c.log 2>&1 &
wait
for t in wb50a wb50b wb50c; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  grep -a "RH-LAND-EAST-ML" $t.log | head -1
  echo "banner diff vs wb49 (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb49.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90
  python3 oceanb.py $V
  python3 landb.py $V | grep -v "Oman\|Namib\|land \|output"
done
touch WB50_DONE
