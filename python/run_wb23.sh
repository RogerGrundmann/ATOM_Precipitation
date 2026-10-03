#!/bin/bash
# 2026-10-03: wb23a / wb23b = wb20b with ATM_MC_GATE_BLEND=3 (parcel carried undiluted, gate threshold 0.01) instead of 1 (user), at
# ATM_RH_OCEAN 0.805 (a, one variable against wb20b) and 0.79 (b, the cliff test). SCREENING: nm 220 from scratch, cli/atm_g3 (-O2), 2 x 8 threads.
# wb20b at 200: 927 mm/a, r .472, sigma 4.12, land/ocean 459/1113, bands 3253/131/210/4.5, land 0-15 1714, Amazon 8.78 mm/d (max cell 74), Congo 2.01;
# tropical land: 14.2 % of cells convect, 2.0 % > 50 mm/d carrying 43 % of the rain; wettest land cell 105 mm/d; ocean east of 128E a uniform 31-32 mm/d.
# Blend mode 1 cliff: wb13a (0.82) 1642 vs wb13b (0.79) 347 mm/a.
# PRE-REGISTERED: every column with 0.01 < |M_u| < 0.1 now rains in proportion to its mass flux, so rain RISES everywhere it was diluted:
# a -- ocean 1113 -> 1500-2200, land 459 -> 700-1300 with the convecting share of tropical land 14 % -> 30-50 % and the share of rain from
# cells > 50 mm/d falling from 43 % to < 25 %; global 1300-1900; sigma BELOW wb20b's relative to the mean (sigma/mean falls);
# b -- the 0.805 -> 0.79 response is smooth: wb23b/wb23a 0.6-0.9 (blend mode 1: 0.21 between 0.82 and 0.79). Deserts stay 0. No NaN;
# RISK: 1/M_u amplification between 0.01 and 0.1 (q_c_u at its 10 g/kg ceiling, single cells > 150 mm/d).
set -u; cd "$(dirname "$0")"; rm -f WB23_DONE
for t in wb23a wb23b; do
  mkdir output_$t || { touch WB23_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB23_DONE; exit 1; }
  sed -e "s#output_wb7/#output_$t/#" -e "s#<nm>600<#<nm>220<#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>220<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_g3 $(md5sum < ../cli/atm_g3 | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 ATM_MC_GATE_BLEND=3 ATM_MC_T_ADD_LAND=2.6 ATM_RH_LAND_QCAP=1"
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN=0.805 ../cli/atm_g3 config_wb23a.xml > wb23a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN=0.79 ../cli/atm_g3 config_wb23b.xml > wb23b.log 2>&1 &
wait
for t in wb23a wb23b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean" $t.log | tail -1
  grep -a "\[MC-GATE\]" $t.log | tail -2 | cut -c1-260
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
  tail -1 output_$t/convergence.csv
  python3 landb.py output_$t/0Ma_smooth_Atm_radial_0_220.vtk
done
touch WB23_DONE
