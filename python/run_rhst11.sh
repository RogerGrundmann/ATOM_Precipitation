#!/bin/bash
# 2026-09-29: STORM -- ATM_RH_STORM=1.1 on the working branch, 600 from scratch, cli/atm_rhst, 8 threads.
# Controls on disk from the same binary and thread count: rhst_ctl (1.0: 35-65 = 69 mm/a), rhst_moist (1.25: 2069).
# PRE-REGISTERED: 35-65 between the two, the response is threshold-like so not linear in the factor; tropics unchanged.
set -u; cd "$(dirname "$0")"; rm -f RHST11_DONE
mkdir output_rhst_11 || { touch RHST11_DONE; exit 1; }
[ -e config_rhst_11.xml ] && { touch RHST11_DONE; exit 1; }
sed "s#output_o2val/#output_rhst_11/#" config_o2val.xml > config_rhst_11.xml
. ./working_branch.env
echo "start $(date +%H:%M)  atm_rhst $(md5sum < ../cli/atm_rhst | cut -c1-8)"
env OMP_NUM_THREADS=8 ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_MFC_DIAG=1 ATM_SR_DIAG=1 ATM_RH_STORM=1.1 \
    ../cli/atm_rhst config_rhst_11.xml > rhst_11.log 2>&1
echo "rhst_11 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' rhst_11.log)  $(date +%H:%M)"
grep -a "by |latitude|" rhst_11.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" rhst_11.log | tail -1 | cut -c1-150
grep -a -i "land .*ocean" rhst_11.log | tail -1 | cut -c1-100; grep -a 'water budget closure' rhst_11.log | tail -1
grep -a 'max v-component' rhst_11.log | tail -1 | cut -c1-80; tail -1 output_rhst_11/convergence.csv
echo "35-65 trajectory: $(grep -a 'by |latitude|' rhst_11.log | grep -v MFC | awk 'NR%6==1{printf "%s ", $12}')"
touch RHST11_DONE
