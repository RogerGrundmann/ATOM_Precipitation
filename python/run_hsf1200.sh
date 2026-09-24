#!/bin/bash
# ATM_HYDRO_SPLIT 600 -> 1200 continuation PAIR from each arm's own iteration-600 checkpoint.
# Binary cli/atm_hsf (c314a4d7), env identical to run_hsf600.sh. 11 threads each, concurrent.
# QUESTION: does the non-hydrostatic p_dyn (0.012 -> 0.191 nd over 600, ~+2e-4/iter, linear)
# PLATEAU? PRE-REGISTERED:
#   1. max|p_dyn| bends over (growth/100 iters falling) and never reaches the 3.0 ceiling
#      -> flip candidate. Stays linear -> investigate the projection, do not flip.
#   2. max|u| tracks hsf1200_ctl; max|v| seam mode not worse than ctl.
#   3. p05 geostrophic residual stays ~0.05; climate null vs ctl at 1200.
#   4. exit 0, zero NaN, zero CEILING BINDING.
set -u; cd "$(dirname "$0")"; rm -f HSF1200_DONE
echo "start $(date +%H:%M)"
( env OMP_NUM_THREADS=11 ATM_RAD_TOPO=0 ATM_UBUD_BALANCE=1 \
      ../cli/atm_hsf config_hsf1200_ctl.xml > hsf1200_ctl.log 2>&1; echo "ctl exit $?" ) &
( env OMP_NUM_THREADS=11 ATM_RAD_TOPO=0 ATM_UBUD_BALANCE=1 ATM_HYDRO_SPLIT=1.0 \
      ../cli/atm_hsf config_hsf1200_on.xml  > hsf1200_on.log  2>&1; echo "on exit $?" ) &
wait
echo "NaN ctl $(grep -c 'NaN/Inf DETECTED' hsf1200_ctl.log)  on $(grep -c 'NaN/Inf DETECTED' hsf1200_on.log)  $(date +%H:%M)"
touch HSF1200_DONE
