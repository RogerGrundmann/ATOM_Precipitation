#!/bin/bash
# THE SNOW WINDOW, CLIMATE ARM (2026-09-27; NOT QUEUED -- launch by hand, only after run_snw40.sh says the deletion
# is large). 600 from scratch, cli/atm_snw, 3 arms x ${NT:-8} threads, ATM_SNOW_DIAG=1 + ATM_CWB_DIAG=1 on all.
#   snw_s0  default (SNOW_WINDOW=0) -- the chain check: must reproduce sgzb_ctl (987 mm/a, r 0.471, sigma 2.37,
#           bands 3388/211/205/24, land/ocean 854/1039) to the thread noise
#   snw_s1  warm side only        snw_s3  both sides
# PRE-REGISTERED against snw_s0: 35-65 deg precipitation rises (the claim under test), 0-15 moves little (warm rain
# dominates there); score on r, sigma, the four bands, land/ocean -- never on the mean alone; the RH_MIN_PTOP fit
# (482) is at stake if the mean moves by more than a few per cent. Exit 0, zero NaN through 155/357/483.
set -u; cd "$(dirname "$0")"; rm -f SNW600_DONE
[ -f SNW40_DONE ] || { echo "run_snw40.sh has not finished"; exit 1; }
for t in s0 s1 s3; do mkdir output_snw_$t || { touch SNW600_DONE; exit 1; }; done
echo "start $(date +%H:%M)"
run(){ ( env OMP_NUM_THREADS=${NT:-8} ATM_SNOW_DIAG=1 ATM_CWB_DIAG=1 ATM_MC_CAP_DIAG=1 ATM_SNOW_WINDOW=$2 ../cli/atm_snw config_snw_$1.xml > snw_$1.log 2>&1
         echo "snw_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' snw_$1.log)  $(date +%H:%M)" ) & }
run s0 0; run s1 1; run s3 3
wait
touch SNW600_DONE
