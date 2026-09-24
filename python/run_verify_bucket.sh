#!/bin/bash
# SHORT CHECKS ONLY (longer runs postponed), 2026-09-24.
# (A) ATM_EVAP_STRIDE_FIX (default 1): the evaporation flux was per ITERATION on a routine called every
#     moist_stride = 2 iterations -> half of E on the ATM_WATER_CLOSURE default. 1 thread, nm = 20.
#       bk_pre  cli/atm_presfix clean      bk_sf0  cli/atm_sfix ATM_EVAP_STRIDE_FIX=0  -- MUST equal bk_pre
#       bk_sfx  cli/atm_sfix clean                                                      -- MUST differ from bk_pre
# (B) ATM_LAND_BUCKET (default 0): bk_bkn cli/atm_bkt unset MUST equal bk_sfx (banner token aside);
#       bk_bk1 cli/atm_bkt ATM_LAND_BUCKET=150 MUST differ.
# (C) 600 -> 620 restart from output_hsf_ctl, ATM_CWB_DIAG=1 ATM_LAND_BUCKET=150 on the default (bkd).
set -u; cd "$(dirname "$0")"; rm -f BK_VERIFY_DONE
. ./verify_lib.sh
for t in pre sf0 sfx bkn bk1; do rm -rf output_bk_$t; mkdir -p output_bk_$t; done
arm bk_pre ../cli/atm_presfix config_bk_pre.xml
arm bk_sf0 ../cli/atm_sfix    config_bk_sf0.xml ATM_EVAP_STRIDE_FIX=0
arm bk_sfx ../cli/atm_sfix    config_bk_sfx.xml
arm bk_bkn ../cli/atm_bkt     config_bk_bkn.xml
arm bk_bk1 ../cli/atm_bkt     config_bk_bk1.xml ATM_LAND_BUCKET=150
( env OMP_NUM_THREADS=3 ATM_CWB_DIAG=1 ATM_LAND_BUCKET=150 ../cli/atm_bkt config_bkd.xml > bkd.log 2>&1 ) &
wait_arms; wait
cmp_dirs    bk_sf0 bk_pre "A1 STRIDE_FIX=0 restores the half-rate branch"
want_differ bk_sfx bk_pre "A2 CONTROL stride fix default vs before -- MUST differ"
cmp_dirs    bk_bkn bk_sfx "B1 LAND_BUCKET unset vs previous build"
want_differ bk_bk1 bk_bkn "B2 CONTROL LAND_BUCKET=150 -- MUST differ"
for p in "sf0 pre" "bkn sfx"; do set -- $p
  diff <(sed "s#output_bk_$1#X#g" output_bk_$1/RUN_CONFIG.txt) <(sed "s#output_bk_$2#X#g" output_bk_$2/RUN_CONFIG.txt) | grep -o "EVAP_STRIDE_FIX=[^ ]*\|LAND_BUCKET=[^ ]*" | tr '\n' ' '; echo; done
echo "== bkd  NaN=$(grep -c 'NaN/Inf DETECTED' bkd.log)"
grep "LAND BUCKET" bkd.log | sed -n '1p;$p' | cut -c1-200
grep "\[CWB\]" bkd.log | grep -E "evaporation|NET|reference" | tail -3 | cut -c1-150
grep -E "model .* NASA|land .* ocean|water budget closure" bkd.log | tail -3 | cut -c1-150
touch BK_VERIFY_DONE
