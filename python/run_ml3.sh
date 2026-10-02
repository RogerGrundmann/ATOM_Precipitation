#!/bin/bash
# 2026-10-02: ATM_MC_ML_PARCEL probe trio -- restart wb7 600 -> 640, working branch, cli/atm_ml (-O2), 3 x 8 threads.
#   ml_ctl  knobs off (= the scheme's saturated cloud-base parcel)
#   ml_A    ML_PARCEL=1, T_ADD 0.15 K, Q_ADD 4.5e-4  (Holtslag-Boville excess at Earth-like ocean fluxes)
#   ml_B    ML_PARCEL=2, T_ADD 0.4 K, dq = sub-grid sigma_q = (1-H_crit) q_sat / sqrt(3) (~1.7 g/kg at 950 hPa)
# PRE-REGISTERED: A -- convection nearly dead everywhere (active < 10 %, P_conv << ctl). B -- buoyant mainly over the
# tropical ocean (ML theta_e + ~4.6 K vs a 3.9 K deficit), weaker over tropical land (6.2 K deficit): ocean<30 P_conv
# similar to or above ctl's 387, land<30 well below ctl's 1784, so the land/ocean ratio falls toward NASA's 0.74.
set -u; cd "$(dirname "$0")"; rm -f ML3_DONE
for t in ml_ctl ml_A ml_B; do
  mkdir output_$t || { touch ML3_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch ML3_DONE; exit 1; }
  ln -s ../output_wb7/atm_restart_0Ma_600.bin output_$t/atm_restart_0Ma_600.bin
  sed -e "s#output_wb7/#output_$t/#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
      -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" -e "s#<nm>600<#<nm>640<#" config_wb7.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_ml $(md5sum < ../cli/atm_ml | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ../cli/atm_ml config_ml_ctl.xml > ml_ctl.log 2>&1 &
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD=0.15 ATM_MC_Q_ADD=4.5e-4 ../cli/atm_ml config_ml_A.xml > ml_A.log 2>&1 &
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_MC_ML_PARCEL=2 ATM_MC_T_ADD=0.4 ../cli/atm_ml config_ml_B.xml > ml_B.log 2>&1 &
wait
for t in ml_ctl ml_A ml_B; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "\[MC-LO\]" $t.log | tail -4 | cut -c1-330
  grep -a "\[MC-PB\]" $t.log | tail -4 | cut -c160-420
  grep -a -i "land .*ocean" $t.log | tail -1 | cut -c1-100
  grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150; grep -a "P_conv mean" $t.log | tail -1
done
touch ML3_DONE
