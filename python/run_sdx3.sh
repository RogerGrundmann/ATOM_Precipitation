#!/bin/bash
# 2026-10-01: third arm of the sdx pair -- ATM_SNOW_DEP_FLUX=1 + ATM_MC_COND_DEBIT=1, working branch (RH_STORM 1.15),
# 600 from scratch, 8 threads, launched ~08:48 while sdx_ctl/sdx_on were at iteration ~150.
# ⚠ BINARY: cli/atm_sdx2 (= cli/atm_cdb, -O2, carries both knobs) -- NOT cli/atm_sdx, which predates
# ATM_MC_COND_DEBIT. The off-branches are byte-identical at -O0 (run_vsdx.sh; run_vcdb.sh pending), but at -O2
# -ffast-math a different binary can differ in last bits, so read this arm against sdx_ctl/sdx_on at the
# thread-noise level, not to the digit.
# PRE-REGISTERED: water budget closes in both schemes; convective latent heating ~+46 W/m2 appears; precipitation
# most likely falls further than sdx_on; bands/r/sigma read, not predicted.
set -u; cd "$(dirname "$0")"; rm -f SDX3_DONE
mkdir output_sdx_both || { touch SDX3_DONE; exit 1; }
[ -e config_sdx_both.xml ] && { touch SDX3_DONE; exit 1; }
sed "s#output_rhst_115/#output_sdx_both/#" config_rhst_115.xml > config_sdx_both.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_sdx2 $(md5sum < ../cli/atm_sdx2 | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_SNOW_DIAG=1 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_SR_DIAG=1 \
    ATM_SNOW_DEP_FLUX=1 ATM_MC_COND_DEBIT=1 ../cli/atm_sdx2 config_sdx_both.xml > sdx_both.log 2>&1
t=both
echo "== sdx_$t  NaN $(grep -c 'NaN/Inf DETECTED' sdx_$t.log)  $(date +%H:%M)"
grep -a "by |latitude|" sdx_$t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" sdx_$t.log | tail -1 | cut -c1-150
grep -a -i "land .*ocean" sdx_$t.log | tail -1 | cut -c1-100; grep -a 'water budget closure' sdx_$t.log | tail -1
grep -a "P_rain mean\|P_snow mean\|P_conv mean" sdx_$t.log | tail -3
grep -a "SNOW DIAG\] GLOBAL" sdx_$t.log | tail -1; grep -a "\[MC-Q\] closure" sdx_$t.log | tail -1
grep -a "as RK4 applies them\|RungeKutta bucket\|\[CWB\] NET\|rest = MC_q" sdx_$t.log | tail -4
grep -a 'precipitable water average' sdx_$t.log | tail -1 | cut -c1-60
grep -a 'max v-component\|max temperature' sdx_$t.log | tail -2 | cut -c1-80; tail -1 output_sdx_$t/convergence.csv
touch SDX3_DONE
