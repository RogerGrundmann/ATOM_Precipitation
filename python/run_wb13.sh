#!/bin/bash
# 2026-10-03: ATM_MC_GATE_BLEND pair, 600 from scratch, cli/atm_gb (-O2), 2 x 8 threads. Base = wb10a's stack (ML_PARCEL 1 T .15 Q 4.5e-4,
# ML_LCL 1, RH_LAND 2) + ATM_MC_GATE_BLEND=1 (updraft only);  wb13a RH_OCEAN 0.82 (cf. wb10a),  wb13b RH_OCEAN 0.79 (cf. wb11).
# wb10a: 1148 mm/a, land/ocean 486/1410, r .407, bands 3883/301/238/4.6.   wb11: 103, 232/52, r -.015, bands 141/1.4/198/4.6.
# [MC-GATE]: ocean<30 gate open in 18.1 % of columns at 0.82, 0.1 % at 0.79; open columns carry 100 % of g_p.
# PRE-REGISTERED: the cliff is gone -- wb13b several hundred mm/a (>> wb11's 103), wb13b/wb13a 0.4-0.9;
# wb13a ABOVE wb10a (the 82 % of ocean<30 columns that were shut now rain a little): 1148 -> 1300-1800, 0-15 band higher still.
set -u; cd "$(dirname "$0")"; rm -f WB13_DONE
for t in wb13a wb13b; do
  mkdir output_$t || { touch WB13_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB13_DONE; exit 1; }
  sed "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_gb $(md5sum < ../cli/atm_gb | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 ATM_MC_GATE_BLEND=1"
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN=0.82 ../cli/atm_gb config_wb13a.xml > wb13a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN=0.79 ../cli/atm_gb config_wb13b.xml > wb13b.log 2>&1 &
wait
for t in wb13a wb13b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean" $t.log | tail -1
  grep -a "\[MC-GATE\]" $t.log | tail -2 | cut -c1-260
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
  tail -1 output_$t/convergence.csv
done
touch WB13_DONE
