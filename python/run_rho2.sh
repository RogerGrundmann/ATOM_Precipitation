#!/bin/bash
# 2026-10-02: ATM_RH_OCEAN probe pair, RE-RUN with the knob on TROPICAL OCEAN only (rho pair was confounded: land got it too) -- 600 from scratch, working branch, cli/atm_rho2 (-O2), 2 x 8 threads concurrent.
#   rho2_82    ATM_RH_OCEAN=0.82 (shipped 0.75), the scheme's cloud-base parcel
#   rho2_82ml  + ATM_MC_ML_PARCEL=2 ATM_MC_T_ADD=0.4 (mixed-layer parcel, sub-grid sigma_q)
# Control: wb7 (same physics, binary differs by print-only probes): 446 mm/a, land/ocean 1000/226, r .225, sigma 1.61,
# ocean<30 P_conv ~387, land<30 ~1784.
# PRE-REGISTERED: ocean<30 BL q +~1.5-2 g/kg, ML theta_e +~4-5 K. rho2_82: ocean rain up (stratiform + convective), land
# unchanged. rho2_82ml: the ML parcel is buoyant over the tropical ocean in a large share of columns (MC-PB "> 0" >> 1 %),
# ocean<30 P_conv > 0 (vs 0 in ml_B); land and 35-65 equal to wb7 (rho2_82) / ml_B-like over land (rho2_82ml). No NaN.
set -u; cd "$(dirname "$0")"; rm -f RHO2_DONE
for t in rho2_82 rho2_82ml; do
  mkdir output_$t || { touch RHO2_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch RHO2_DONE; exit 1; }
  sed "s#output_wb7/#output_$t/#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_rho2 $(md5sum < ../cli/atm_rho2 | cut -c1-8)"
D="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1"
env OMP_NUM_THREADS=8 $D ATM_RH_OCEAN=0.82 ../cli/atm_rho2 config_rho2_82.xml > rho2_82.log 2>&1 &
env OMP_NUM_THREADS=8 $D ATM_RH_OCEAN=0.82 ATM_MC_ML_PARCEL=2 ATM_MC_T_ADD=0.4 ../cli/atm_rho2 config_rho2_82ml.xml > rho2_82ml.log 2>&1 &
wait
for t in rho2_82 rho2_82ml; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a -i "land .*ocean" $t.log | tail -1 | cut -c1-100; grep -a 'water budget closure' $t.log | tail -1
  grep -a "P_conv mean" $t.log | tail -1
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
  grep -a "\[MC-PB\] ocean<30\|\[MC-PB\] land<30" $t.log | tail -2
  grep -a "\[BL-Q\]" $t.log | tail -4 | cut -c1-260
  grep -a 'max v-component\|max w-component\|mean temperature' $t.log | tail -3 | cut -c1-80; tail -1 output_$t/convergence.csv
done
touch RHO2_DONE
