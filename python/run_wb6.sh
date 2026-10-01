#!/bin/bash
# 2026-10-01: FULL working branch + ATM_MC_QC_DETRAIN (all knobs at d851438 + QC_DETRAIN), 600 from scratch, cli/atm_qc
# (-O2, f90c462 source), 8 threads. Reference: wb5 = this minus ATM_MC_QC_DETRAIN (cli/atm_tv3): 476.0 mm/a, r +0.363,
# sigma 1.43, bands 1281/260/216/4.6, land/ocean 632/414, P_conv 394.
# PRE-REGISTERED: P_conv falls ~2/3 (qc probe 311 -> 101 over 40 iters), total precip falls by about that loss, to
# ~200-250 mm/a; the tropics 0-15 carry it; extratropical bands unchanged; MC_q = -P_conv; no NaN.
set -u; cd "$(dirname "$0")"; rm -f WB6_DONE
mkdir output_wb6 || { touch WB6_DONE; exit 1; }
[ -e config_wb6.xml ] && { touch WB6_DONE; exit 1; }
sed "s#output_rhst_115/#output_wb6/#" config_rhst_115.xml > config_wb6.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_qc $(md5sum < ../cli/atm_qc | cut -c1-8)  SNOW_DEP_FLUX=$ATM_SNOW_DEP_FLUX MC_COND_DEBIT=$ATM_MC_COND_DEBIT MC_Q_NDIM=$ATM_MC_Q_NDIM Q_DIFF_FLUX=$ATM_Q_DIFF_FLUX MC_T_COEFF=$ATM_MC_T_COEFF MC_UV_DETRAIN=$ATM_MC_UV_DETRAIN MC_VEL_COEFF=$ATM_MC_VEL_COEFF MC_QC_DETRAIN=$ATM_MC_QC_DETRAIN RH_STORM=$ATM_RH_STORM"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_SR_DIAG=1 ATM_SNOW_DIAG=1 \
    ../cli/atm_qc config_wb6.xml > wb6.log 2>&1
t=wb6
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
grep -a -i "land .*ocean" $t.log | tail -1 | cut -c1-100; grep -a 'water budget closure' $t.log | tail -1
grep -a "P_rain mean\|P_snow mean\|P_conv mean" $t.log | tail -3
grep -a "SNOW DIAG\] GLOBAL" $t.log | tail -1; grep -a "\[MC-Q\] closure" $t.log | tail -1
grep -a "\[CWB\] evaporation\|\[CWB\] NET\|RungeKutta bucket\|rest = MC_q\|CWB-TD\]   advection + diff" $t.log | tail -5
grep -a 'precipitable water average' $t.log | tail -1 | cut -c1-60
grep -a 'max v-component\|max w-component\|max w_u \|mean temperature\|max temperature' $t.log | tail -5 | cut -c1-80; grep -a '\[MC-TV\]\|\[MC-QC\]\|\[MC DIAG\] column budget' $t.log | tail -5; tail -1 output_$t/convergence.csv
touch WB6_DONE
