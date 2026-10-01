#!/bin/bash
# 2026-10-01: FULL working branch at 105e391 (closure + convection stack + RH_STORM 1.15 + SNOW_DEP_FLUX + MC_COND_DEBIT +
# MC_Q_NDIM + Q_DIFF_FLUX), 600 from scratch, cli/atm_qdf (-O2, all four knobs), 8 threads.
# Reference: wb3 = this minus ATM_Q_DIFF_FLUX (cli/atm_qnd): 475.7 mm/a, r +0.363, sigma 1.43, bands 1279/260/216/4.6,
# land/ocean 633/414. PRE-REGISTERED: CWB transport/diffusion ~ -14 (the advection); ~220 mm/a less water made aloft,
# so precipitation most likely FALLS and the column moistens less; bands/r/sigma read, not predicted.
set -u; cd "$(dirname "$0")"; rm -f WB4_DONE
mkdir output_wb4 || { touch WB4_DONE; exit 1; }
[ -e config_wb4.xml ] && { touch WB4_DONE; exit 1; }
sed "s#output_rhst_115/#output_wb4/#" config_rhst_115.xml > config_wb4.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_qdf $(md5sum < ../cli/atm_qdf | cut -c1-8)  SNOW_DEP_FLUX=$ATM_SNOW_DEP_FLUX MC_COND_DEBIT=$ATM_MC_COND_DEBIT MC_Q_NDIM=$ATM_MC_Q_NDIM Q_DIFF_FLUX=$ATM_Q_DIFF_FLUX RH_STORM=$ATM_RH_STORM"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_SR_DIAG=1 ATM_SNOW_DIAG=1 \
    ../cli/atm_qdf config_wb4.xml > wb4.log 2>&1
t=wb4
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
grep -a -i "land .*ocean" $t.log | tail -1 | cut -c1-100; grep -a 'water budget closure' $t.log | tail -1
grep -a "P_rain mean\|P_snow mean\|P_conv mean" $t.log | tail -3
grep -a "SNOW DIAG\] GLOBAL" $t.log | tail -1; grep -a "\[MC-Q\] closure" $t.log | tail -1
grep -a "\[CWB\] evaporation\|\[CWB\] NET\|RungeKutta bucket\|rest = MC_q\|CWB-TD\]   advection + diff" $t.log | tail -5
grep -a 'precipitable water average' $t.log | tail -1 | cut -c1-60
grep -a 'max v-component' $t.log | tail -1 | cut -c1-80; tail -1 output_$t/convergence.csv
touch WB4_DONE
