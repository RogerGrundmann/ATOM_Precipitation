#!/bin/bash
# 2026-10-06: wb48a / b / c = wb47b (working branch + ATM_RH_LAND_EAST=1.6 + ATM_RH_LAND_EAST_MAX) with a LOOSER cap: +0.02 / +0.03 / +0.05 (wb47b: +0.010).
# User 10-06: "adopt wb47b and run the 600, screen a looser cap first" -- this is the screen (the trade curve cap -> E Australia vs land 15-35).
# SCREENING AT nm 60 from scratch, cli/atm_emx (-O2, md5 3fef9110, same binary as wb47), 3 x 8 threads, ATM_VTK_STRIDE=2. One variable against wb47b.
# Control = wb47b (60): 971.8 mm/a (-0.7 %), r .609, sigma 1.24, land 0-15 / 15-35 / 35-65 1664 / 217 / 668, E Austral 2.73 (NASA 1.86), Madagascar 3.16 (4.31),
# S China 0.81 (3.87), SE US 1.19 (3.89), S Brazil 0.85 (4.55), Arabia 0.21, wettest cell 10.8 mm/d.
# PRE-REGISTERED (at 60; the rain-vs-RH response is a cliff -- lowland 15-35: RH 0.80-0.84 2.46, 0.84-0.88 6.68, > 0.88 10.5 mm/d -- so wide ranges):
#   a +0.02: E Austral 3.3-4.5, land 15-35 250-320, wettest cell 11-14;  b +0.03: E Austral 4-6, land 15-35 300-420, wettest cell 12-16;
#   c +0.05: E Austral 5.5-9, land 15-35 400-700, wettest cell > 16, r below .60.  S China / SE US rise only where their cells sat ON the cap.
# QUESTION: is there a cap with land 15-35 > 320 and E Australia < 2x NASA (3.7)? If not, wb47b stands.
# RESULT (08:11-08:16, NaN 0; land 15-35 / E Austral / Madagascar / S China / SE US / S Brazil / global / r / sigma / wettest cell):
#   wb47b +0.010: 217 / 2.73 / 3.16 / 0.81 / 1.19 / 0.85 / -0.7 % / .609 / 1.24 / 10.8
#   a +0.02:      276 / 3.90 / 3.99 / 1.20 / 1.56 / 1.37 / +0.1 % / .610 / 1.25 / 12.6
#   b +0.03:      342 / 5.14 / 5.16 / 1.64 / 1.93 / 2.05 / +0.8 % / .609 / 1.26 / 13.7
#   c +0.05:      477 / 8.09 / 7.82 / 2.64 / 2.61 / 3.68 / +2.2 % / .602 / 1.29 / 18.1
# ANSWER: no. Land 15-35 passes 320 only from +0.03, where E Australia is already 2.8x NASA; every east side rises together (one humidity, one cliff).
# r and sigma are flat up to +0.03. wb47b stands and goes to 600 (run_wb49.sh); +0.02 is the alternative if S China / SE US count more than E Australia.
set -u; cd "$(dirname "$0")"; rm -f WB48_DONE
for t in wb48a wb48b wb48c; do
  mkdir output_$t || { touch WB48_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB48_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>60<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>60<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_emx $(md5sum < ../cli/atm_emx | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_VTK_STRIDE=2 ATM_RH_LAND_EAST=1.6"
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST_MAX=0.02 ../cli/atm_emx config_wb48a.xml > wb48a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST_MAX=0.03 ../cli/atm_emx config_wb48b.xml > wb48b.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_LAND_EAST_MAX=0.05 ../cli/atm_emx config_wb48c.xml > wb48c.log 2>&1 &
wait
for t in wb48a wb48b wb48c; do
  V=output_$t/0Ma_smooth_Atm_radial_0_60.vtk
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)  slices: $(ls output_$t | grep -c radial_0)"
  echo "banner diff vs wb47b (knob tokens only):"; diff <(grep -a "RUN CONFIG" wb47b.log | tr ' ' '\n' | grep "=" | sort -u) <(grep -a "RUN CONFIG" $t.log | tr ' ' '\n' | grep "=" | sort -u) | grep "^[<>]" | grep -v "nm=\|output" | head -8
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90
  python3 oceanb.py $V
  python3 landb.py $V | grep -v "Oman\|Namib\|land \|output"
done
touch WB48_DONE
