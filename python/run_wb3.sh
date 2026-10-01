#!/bin/bash
# 2026-10-01: the FULL working branch after the water-budget fixes (working_branch.env at e0bcbba: closure stack +
# convection stack + RH_STORM 1.15 + SNOW_DEP_FLUX + MC_COND_DEBIT + MC_Q_NDIM), 600 from scratch, cli/atm_qnd (-O2,
# carries all three knobs), 8 threads. Runs alongside rs2_100/rs2_125 (3 x 8 threads).
# Reference: sdx_both = this minus ATM_MC_Q_NDIM (other binary, cli/atm_sdx2): 412.1 mm/a, r +0.316, sigma 1.24,
# bands 1011/272/216/4.6, land/ocean 645/320.
# PRE-REGISTERED: climate null against sdx_both (thread noise); CWB MC_q row = -P_conv; NET ~ E - P in the surface band.
set -u; cd "$(dirname "$0")"; rm -f WB3_DONE
mkdir output_wb3 || { touch WB3_DONE; exit 1; }
[ -e config_wb3.xml ] && { touch WB3_DONE; exit 1; }
sed "s#output_rhst_115/#output_wb3/#" config_rhst_115.xml > config_wb3.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_qnd $(md5sum < ../cli/atm_qnd | cut -c1-8)  SNOW_DEP_FLUX=$ATM_SNOW_DEP_FLUX MC_COND_DEBIT=$ATM_MC_COND_DEBIT MC_Q_NDIM=$ATM_MC_Q_NDIM RH_STORM=$ATM_RH_STORM"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_SR_DIAG=1 ATM_SNOW_DIAG=1 \
    ../cli/atm_qnd config_wb3.xml > wb3.log 2>&1
t=wb3
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
grep -a -i "land .*ocean" $t.log | tail -1 | cut -c1-100; grep -a 'water budget closure' $t.log | tail -1
grep -a "P_rain mean\|P_snow mean\|P_conv mean" $t.log | tail -3
grep -a "SNOW DIAG\] GLOBAL" $t.log | tail -1; grep -a "\[MC-Q\] closure" $t.log | tail -1
grep -a "\[CWB\] evaporation\|\[CWB\] NET\|RungeKutta bucket\|rest = MC_q" $t.log | tail -4
grep -a 'precipitable water average' $t.log | tail -1 | cut -c1-60
grep -a 'max v-component' $t.log | tail -1 | cut -c1-80; tail -1 output_$t/convergence.csv
touch WB3_DONE
