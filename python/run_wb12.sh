#!/bin/bash
# 2026-10-02: ATM_MC_DEPTH_RAMP pair, 600 from scratch, cli/atm_dr (-O2), 2 x 8 threads. Base = wb10a's stack (ML_PARCEL 1 T .15 Q 4.5e-4,
# ML_LCL 1, RH_LAND 2) + ATM_MC_DEPTH_RAMP=1;  wb12a RH_OCEAN 0.82 (cf. wb10a),  wb12b RH_OCEAN 0.79 (cf. wb11).
# wb10a: 1148 mm/a, land/ocean 486/1410, r .407, bands 3883/301/238/4.6.   wb11: 103, 232/52, r -.015, bands 141/1.4/198/4.6.
# PRE-REGISTERED: the cliff is gone -- wb12b well above wb11 (several hundred mm/a, ocean<15 g_p >> 18); wb12a close to wb10a
# (shallow columns now rain a little, 200-300 hPa columns a little less); the 0.82 -> 0.79 response is smooth (wb12b/wb12a ~0.5-0.9).
set -u; cd "$(dirname "$0")"; rm -f WB12_DONE
for t in wb12a wb12b; do
  mkdir output_$t || { touch WB12_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch WB12_DONE; exit 1; }
  sed "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_dr $(md5sum < ../cli/atm_dr | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ATM_MC_ML_LCL=1 ATM_RH_LAND=2 ATM_MC_DEPTH_RAMP=1"
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN=0.82 ../cli/atm_dr config_wb12a.xml > wb12a.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_RH_OCEAN=0.79 ../cli/atm_dr config_wb12b.xml > wb12b.log 2>&1 &
wait
for t in wb12a wb12b; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean" $t.log | tail -1
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-RG\]" $t.log | tail -6 | cut -c1-260
  tail -1 output_$t/convergence.csv
done
touch WB12_DONE
