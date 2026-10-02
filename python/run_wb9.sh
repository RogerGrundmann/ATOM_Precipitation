#!/bin/bash
# 2026-10-02: wb9 / wb9n, 600 from scratch, cli/atm_rhl (-O2), 2 x 8 threads. Working branch + ATM_RH_OCEAN=0.82 + ATM_MC_ML_PARCEL=1
# (fixed excess, Holtslag-Boville Earth-like: T_ADD 0.15 K, Q_ADD 4.5e-4) + ATM_MC_ML_LCL=1;  wb9 adds ATM_RH_LAND=1 (land initial RH from
# subtropical descent x continentality, paleo-safe). wb9n isolates RH_LAND; wb9n vs wb8 isolates the excess mode (1 vs 2).
# wb8 (ML_PARCEL 2, T 0.4): 694 mm/a, land/ocean 798/652, r .158, sigma 4.35; land hotspots Horn/Arabia ~130 mm/d.
# PRE-REGISTERED: wb9n -- the elevated-desert hotspots weaken (no pressure-dependent sigma_q) but deserts still rain (humid
# BL); ocean convection weaker than wb8 (smaller excess). wb9 -- Sahara/Arabia/Australia rain falls toward NASA (<1 mm/d),
# Amazon/Congo unchanged vs wb9n, share of land rain from cells > 3x NASA falls well below wb8's ~84 %, sigma falls, r rises.
set -u; cd "$(dirname "$0")"; rm -f WB9_DONE
for t in wb9 wb9n; do
  mkdir output_$t || { touch WB9_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB9_DONE; exit 1; }
  sed "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_rhl $(md5sum < ../cli/atm_rhl | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_RH_OCEAN=0.82 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1"
env OMP_NUM_THREADS=8 $K ATM_RH_LAND=1 ../cli/atm_rhl config_wb9.xml  > wb9.log  2>&1 &
env OMP_NUM_THREADS=8 $K              ../cli/atm_rhl config_wb9n.xml > wb9n.log 2>&1 &
wait
for t in wb9 wb9n; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "RH-LAND" $t.log | head -1
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean" $t.log | tail -1
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  tail -1 output_$t/convergence.csv
done
touch WB9_DONE
