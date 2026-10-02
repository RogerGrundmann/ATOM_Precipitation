#!/bin/bash
# 2026-10-02: ATM_MC_ML_LCL probe pair -- restart rho2_82ml 600 -> 640, working branch + RH_OCEAN 0.82 + ML_PARCEL 2 + T_ADD 0.4,
# cli/atm_lcl (-O2), 2 x 8 threads.   lcl_0: ML_LCL=0 (= rho2_82ml)   lcl_1: ML_LCL=1 (base = the ML parcel's LCL)
# PRE-REGISTERED: lcl_1 ocean<30 base 1722 m -> ~500-800 m, q_c_u(base) 2.4 -> < 0.5 g/kg, P_conv max << 1.2e5 mm/a,
# sigma falls well below 9.4; ocean<30 P_conv stays > 0 (ML parcel buoyant in ~53 % of columns); land<30 P_conv falls
# (land ML parcel buoyant in only ~5 %). No NaN.
set -u; cd "$(dirname "$0")"; rm -f LCL_DONE
for t in lcl_0 lcl_1; do
  mkdir output_$t || { touch LCL_DONE; exit 1; }
  [ -e config_$t.xml ] && { touch LCL_DONE; exit 1; }
  ln -s ../output_rho2_82ml/atm_restart_0Ma_600.bin output_$t/atm_restart_0Ma_600.bin
  sed -e "s#output_rho2_82ml/#output_$t/#" -e "s#<checkpoint_save_iter>600<#<checkpoint_save_iter>-1<#" \
      -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" -e "s#<nm>600<#<nm>640<#" config_rho2_82ml.xml > config_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_lcl $(md5sum < ../cli/atm_lcl | cut -c1-8)"
K="ATM_MC_DIAG=1 ATM_RH_OCEAN=0.82 ATM_MC_ML_PARCEL=2 ATM_MC_T_ADD=0.4"
env OMP_NUM_THREADS=8 $K ATM_MC_ML_LCL=0 ../cli/atm_lcl config_lcl_0.xml > lcl_0.log 2>&1 &
env OMP_NUM_THREADS=8 $K ATM_MC_ML_LCL=1 ../cli/atm_lcl config_lcl_1.xml > lcl_1.log 2>&1 &
wait
for t in lcl_0 lcl_1; do
  echo "== $t  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
  grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
  grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a "P_conv mean" $t.log | tail -1
  grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
done
touch LCL_DONE
