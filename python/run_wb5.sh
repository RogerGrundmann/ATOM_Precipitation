#!/bin/bash
# 2026-10-01: FULL working branch + the three MC-TV knobs (ATM_MC_T_COEFF, ATM_MC_UV_DETRAIN, ATM_MC_VEL_COEFF), 600 from
# scratch, cli/atm_tv3 (-O2, carries all knobs), 8 threads. Reference: wb4 = this minus the three MC-TV knobs (cli/atm_qdf):
# 476.9 mm/a, r +0.364, sigma 1.43, bands 1285/260/216/4.6, land/ocean 631/416.
# PRE-REGISTERED: mean T unchanged (teq relaxation eats the +24 W/m2 heating); precip and bands null to ~1 %;
# the zonal wind where convection acts weakly smoothed (drag/cor ~0.09) -- max winds within a few %; no NaN.
set -u; cd "$(dirname "$0")"; rm -f WB5_DONE
mkdir output_wb5 || { touch WB5_DONE; exit 1; }
[ -e config_wb5.xml ] && { touch WB5_DONE; exit 1; }
sed "s#output_rhst_115/#output_wb5/#" config_rhst_115.xml > config_wb5.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_tv3 $(md5sum < ../cli/atm_tv3 | cut -c1-8)  SNOW_DEP_FLUX=$ATM_SNOW_DEP_FLUX MC_COND_DEBIT=$ATM_MC_COND_DEBIT MC_Q_NDIM=$ATM_MC_Q_NDIM Q_DIFF_FLUX=$ATM_Q_DIFF_FLUX MC_T_COEFF=$ATM_MC_T_COEFF MC_UV_DETRAIN=$ATM_MC_UV_DETRAIN MC_VEL_COEFF=$ATM_MC_VEL_COEFF RH_STORM=$ATM_RH_STORM"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_SR_DIAG=1 ATM_SNOW_DIAG=1 \
    ../cli/atm_tv3 config_wb5.xml > wb5.log 2>&1
t=wb5
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
grep -a -i "land .*ocean" $t.log | tail -1 | cut -c1-100; grep -a 'water budget closure' $t.log | tail -1
grep -a "P_rain mean\|P_snow mean\|P_conv mean" $t.log | tail -3
grep -a "SNOW DIAG\] GLOBAL" $t.log | tail -1; grep -a "\[MC-Q\] closure" $t.log | tail -1
grep -a "\[CWB\] evaporation\|\[CWB\] NET\|RungeKutta bucket\|rest = MC_q\|CWB-TD\]   advection + diff" $t.log | tail -5
grep -a 'precipitable water average' $t.log | tail -1 | cut -c1-60
grep -a 'max v-component\|max w-component\|max w_u \|mean temperature\|max temperature' $t.log | tail -5 | cut -c1-80; grep -a '\[MC-TV\]' $t.log | tail -3; tail -1 output_$t/convergence.csv
touch WB5_DONE
