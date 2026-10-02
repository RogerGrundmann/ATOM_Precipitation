#!/bin/bash
# 2026-10-02: ATM_RH_OCEAN probe pair -- 600 from scratch, working branch, cli/atm_rho (-O2), 2 x 8 threads concurrent.
#   rho_82    ATM_RH_OCEAN=0.82 (shipped 0.75), the scheme's cloud-base parcel
#   rho_82ml  + ATM_MC_ML_PARCEL=2 ATM_MC_T_ADD=0.4 (mixed-layer parcel, sub-grid sigma_q)
# Control: wb7 (same physics, binary differs by print-only probes): 446 mm/a, land/ocean 1000/226, r .225, sigma 1.61,
# ocean<30 P_conv ~387, land<30 ~1784.
# PRE-REGISTERED: ocean<30 BL q +~1.5-2 g/kg, ML theta_e +~4-5 K. rho_82: ocean rain up (stratiform + convective), land
# unchanged. rho_82ml: the ML parcel is buoyant over the tropical ocean in a large share of columns (MC-PB "> 0" >> 1 %),
# ocean<30 P_conv > 0 (vs 0 in ml_B); land<30 ~ ml_B (~2000). No NaN.
set -u; cd "$(dirname "$0")"; rm -f RHO_DONE
for t in rho_82 rho_82ml; do
  mkdir output_$t || { touch RHO_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch RHO_DONE; exit 1; }
  sed "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_rho $(md5sum < ../cli/atm_rho | cut -c1-8)"
D="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1"
env OMP_NUM_THREADS=8 $D ATM_RH_OCEAN=0.82 ../cli/atm_rho config_rho_82.xml > rho_82.log 2>&1 &
env OMP_NUM_THREADS=8 $D ATM_RH_OCEAN=0.82 ATM_MC_ML_PARCEL=2 ATM_MC_T_ADD=0.4 ../cli/atm_rho config_rho_82ml.xml > rho_82ml.log 2>&1 &
wait
for t in rho_82 rho_82ml; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a -i "land .*ocean" $t.log | tail -1 | cut -c1-100; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean" $t.log | tail -1
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-PB\] ocean<30\|\[MC-PB\] land<30" $t.log | tail -2
  grep -a "\[BL-Q\]" $t.log | tail -4 | cut -c1-260
  grep -a 'max v-component\|max w-component\|mean temperature' $t.log | tail -3 | cut -c1-80; tail -1 output_$t/convergence.csv
done
touch RHO_DONE
