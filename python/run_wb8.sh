#!/bin/bash
# 2026-10-02: wb8 = working branch + ATM_RH_OCEAN=0.82 + ATM_MC_ML_PARCEL=2 + ATM_MC_T_ADD=0.4 + ATM_MC_ML_LCL=1, 600 from scratch,
# cli/atm_lcl (-O2), 8 threads. References: wb7 446 mm/a, land/ocean 1000/226, r .225, sigma 1.61; rho2_82ml (no LCL) 1037,
# 1037/1038, r .179, sigma 9.37. lcl_1 probe (40 iters from rho2_82ml): 689, land/ocean 795/646, r .155, sigma 4.32.
# PRE-REGISTERED: from scratch the convective base sits at the ML parcel's LCL (~0.8 km ocean<30), no dumped base condensate
# (q_c_u(base) < 0.3 g/kg), sigma between 1.6 and 4.3, land ~800, ocean 500-900, global 600-800. No NaN.
set -u; cd "$(dirname "$0")"; rm -f WB8_DONE
mkdir output_wb8 || { touch WB8_DONE; exit 1; }
[ -e config_wb8.xml ] && { touch WB8_DONE; exit 1; }
sed "s#output_wb7/#output_wb8/#" config_wb7.xml > config_wb8.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_lcl $(md5sum < ../cli/atm_lcl | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_RH_OCEAN=0.82 ATM_MC_ML_PARCEL=2 ATM_MC_T_ADD=0.4 ATM_MC_ML_LCL=1 \
    ../cli/atm_lcl config_wb8.xml > wb8.log 2>&1
t=wb8
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
grep -a 'land .*ocean .*(model / NASA)' $t.log | tail -1 | cut -c1-90; grep -a 'water budget closure' $t.log | tail -1
grep -a "P_conv mean" $t.log | tail -1
grep -a "\[MC-LO\] ocean<30\|\[MC-LO\] land<30" $t.log | tail -2 | cut -c1-330
grep -a "\[MC-PB\] ocean<30\|\[MC-PB\] land<30" $t.log | tail -2
grep -a 'max v-component\|max w-component\|mean temperature' $t.log | tail -3 | cut -c1-80; tail -1 output_$t/convergence.csv
touch WB8_DONE
