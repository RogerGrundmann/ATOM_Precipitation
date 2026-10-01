#!/bin/bash
# 2026-10-01: FULL working branch + the convective parcel condensing (ATM_MC_SGZ=1, ATM_MC_ENTR=1e-4) + ATM_MC_ED_ABOVE_BASE=1,
# 600 from scratch, cli/atm_eab (-O2, 8c8d7bd source), 8 threads. Reference: wb6 = this minus those three (cli/atm_qc):
# 228.3 mm/a, r +0.326, sigma 0.71, bands 531/71/216/4.6, land/ocean 213/234, P_conv 146.
# PRE-REGISTERED: P_conv rises well above wb6 (pcd2 40-iter: 101 -> 330), total ~400-450 mm/a, carried by 0-15 and 15-35;
# 35-65/65-90 unchanged; MC_q = -P_conv; no NaN. Unknown: whether the extra sub-cloud moistening feeds back.
set -u; cd "$(dirname "$0")"; rm -f WB7_DONE
mkdir output_wb7 || { touch WB7_DONE; exit 1; }
[ -e config_wb7.xml ] && { touch WB7_DONE; exit 1; }
sed "s#output_rhst_115/#output_wb7/#" config_rhst_115.xml > config_wb7.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_eab $(md5sum < ../cli/atm_eab | cut -c1-8)  SNOW_DEP_FLUX=$ATM_SNOW_DEP_FLUX MC_COND_DEBIT=$ATM_MC_COND_DEBIT MC_Q_NDIM=$ATM_MC_Q_NDIM Q_DIFF_FLUX=$ATM_Q_DIFF_FLUX MC_T_COEFF=$ATM_MC_T_COEFF MC_UV_DETRAIN=$ATM_MC_UV_DETRAIN MC_VEL_COEFF=$ATM_MC_VEL_COEFF MC_QC_DETRAIN=$ATM_MC_QC_DETRAIN SGZ=$ATM_MC_SGZ ENTR=$ATM_MC_ENTR ED_ABOVE_BASE=$ATM_MC_ED_ABOVE_BASE RH_STORM=$ATM_RH_STORM"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_SR_DIAG=1 ATM_SNOW_DIAG=1 \
    ../cli/atm_eab config_wb7.xml > wb7.log 2>&1
t=wb7
echo "== $t  exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $t.log)  $(date +%H:%M)"
grep -a "by |latitude|" $t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" $t.log | tail -1 | cut -c1-150
grep -a -i "land .*ocean" $t.log | tail -1 | cut -c1-100; grep -a 'water budget closure' $t.log | tail -1
grep -a "P_rain mean\|P_snow mean\|P_conv mean" $t.log | tail -3
grep -a "SNOW DIAG\] GLOBAL" $t.log | tail -1; grep -a "\[MC-Q\] closure" $t.log | tail -1
grep -a "\[CWB\] evaporation\|\[CWB\] NET\|RungeKutta bucket\|rest = MC_q\|CWB-TD\]   advection + diff" $t.log | tail -5
grep -a 'precipitable water average' $t.log | tail -1 | cut -c1-60
grep -a 'max v-component\|max w-component\|max w_u \|mean temperature\|max temperature' $t.log | tail -5 | cut -c1-80; grep -a '\[MC-TV\]\|\[MC-QC\]\|\[MC DIAG\] column budget' $t.log | tail -5; tail -1 output_$t/convergence.csv
touch WB7_DONE
