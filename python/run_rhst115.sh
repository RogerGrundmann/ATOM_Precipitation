#!/bin/bash
# 2026-09-30: STORM -- ATM_RH_STORM=1.15 (the adopted working-branch value) on the working branch, 600 from scratch, cli/atm_rhst, 8 threads.
# Controls on disk from the same binary and thread count: rhst_ctl (1.0: 35-65 = 69 mm/a), rhst_moist (1.25: 2069), rhst_11 (1.1: 313).
# PRE-REGISTERED: 35-65 between 313 and 2069, the response is threshold-like so not linear in the factor; tropics unchanged.
set -u; cd "$(dirname "$0")"; rm -f RHST115_DONE
mkdir output_rhst_115 || { touch RHST115_DONE; exit 1; }
[ -e config_rhst_115.xml ] && { touch RHST115_DONE; exit 1; }
sed "s#output_o2val/#output_rhst_115/#" config_o2val.xml > config_rhst_115.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_rhst $(md5sum < ../cli/atm_rhst | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_MFC_DIAG=1 ATM_SR_DIAG=1 ATM_RH_STORM=1.15 \
    ../cli/atm_rhst config_rhst_115.xml > rhst_115.log 2>&1
echo "rhst_115 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' rhst_115.log)  $(date +%H:%M)"
grep -a "by |latitude|" rhst_115.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" rhst_115.log | tail -1 | cut -c1-150
grep -a -i "land .*ocean" rhst_115.log | tail -1 | cut -c1-100; grep -a 'water budget closure' rhst_115.log | tail -1
grep -a 'max v-component' rhst_115.log | tail -1 | cut -c1-80; tail -1 output_rhst_115/convergence.csv
echo "35-65 trajectory: $(grep -a 'by |latitude|' rhst_115.log | grep -v MFC | awk 'NR%6==1{printf "%s ", $12}')"
touch RHST115_DONE
