#!/bin/bash
# 2026-10-01: STORM re-judged on the corrected working branch (ATM_SNOW_DEP_FLUX=1 + ATM_MC_COND_DEBIT=1, now in
# working_branch.env). ATM_RH_STORM 1.0 and 1.25, 600 from scratch, cli/atm_sdx2 (-O2; the same binary as sdx_both,
# which IS the 1.15 point), 2 x 8 threads concurrent.
# On the leaking branch: 1.0 -> 35-65 69, 1.1 -> 313, 1.15 -> 632, 1.25 -> 2069 (run_rhst*.sh). With the leak
# removed 1.15 gives 216. PRE-REGISTERED: monotone in the factor again, every point well below its leaking-branch
# value; 1.25 no longer breaks the shape the way 2069 did; whether ANY value reaches 981 without overshooting land
# is the question, not predicted.
set -u; cd "$(dirname "$0")"; rm -f RS2_DONE
for t in 100 125; do
  mkdir output_rs2_$t || { touch RS2_DONE; exit 1; }
  [ -e config_rs2_$t.xml ] && { touch RS2_DONE; exit 1; }
  sed "s#output_rhst_115/#output_rs2_$t/#" config_rhst_115.xml > config_rs2_$t.xml
done
. ./working_branch.env
echo "start $(date +%H:%M)  atm_sdx2 $(md5sum < ../cli/atm_sdx2 | cut -c1-8)  SNOW_DEP_FLUX=$ATM_SNOW_DEP_FLUX MC_COND_DEBIT=$ATM_MC_COND_DEBIT"
D="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_SR_DIAG=1 ATM_SNOW_DIAG=1"
env OMP_NUM_THREADS=8 $D ATM_RH_STORM=1.0  ../cli/atm_sdx2 config_rs2_100.xml > rs2_100.log 2>&1 &
env OMP_NUM_THREADS=8 $D ATM_RH_STORM=1.25 ../cli/atm_sdx2 config_rs2_125.xml > rs2_125.log 2>&1 &
wait
for t in 100 125; do
  echo "== rs2_$t  NaN $(grep -c 'NaN/Inf DETECTED' rs2_$t.log)  $(date +%H:%M)  $(grep -ao 'RH_STORM=[^ ]*' output_rs2_$t/RUN_CONFIG.txt | head -1)"
  grep -a "by |latitude|" rs2_$t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" rs2_$t.log | tail -1 | cut -c1-150
  grep -a -i "land .*ocean" rs2_$t.log | tail -1 | cut -c1-100; grep -a 'water budget closure' rs2_$t.log | tail -1
  grep -a "P_rain mean\|P_snow mean\|P_conv mean" rs2_$t.log | tail -3
  grep -a "\[CWB\] NET" rs2_$t.log | tail -1; grep -a 'precipitable water average' rs2_$t.log | tail -1 | cut -c1-60
  tail -1 output_rs2_$t/convergence.csv
  echo "35-65 trajectory: $(grep -a 'by |latitude|' rs2_$t.log | grep -v MFC | awk 'NR%6==1{printf "%s ", $12}')"
done
touch RS2_DONE
