#!/bin/bash
# 2026-09-29: STORM -- is the closure branch's starved storm track the SaturationAdjustment latent heat that sync
# mode 2 keeps? On the shipped branch RK4 discards SatAdj's direct writes (only ~0.02-0.8 % survive via S_c_c);
# ATM_WATER_CLOSURE forces ATM_RK_SCALAR_SYNC=2, which keeps both the condensate AND its heat. qh_on's 35-65 column
# was +0.65 K warmer than the closure-off arms. rksq_on (100 iters, sync 1) was the first probe and is unscored.
# Setup = qh_on (closure + moisture filter off + vertical T filter off) ON THE CONVERGED RAIN BRANCH
# (ATM_PRECIP_UPWIND=1 ATM_SNOW_WINDOW=2), 600 from scratch, cli/atm_sdr (same binary as the RAIN_AREA sweep,
# pre-SURF_DRAG-flip default 0.0 -- measured climate-null), 2 x 6 threads = the sweep's per-arm thread count.
#   stm_2  closure, sync 2 (moisture + t)  -- the closure branch as it ships
#   stm_1  closure, ATM_RK_SCALAR_SYNC=1 (moisture only: condensate kept, its latent heat discarded)
# control on the shipped branch: ra_10 (run_ra.sh: RAIN_AREA 0.10, same UPWIND + SNOW_WINDOW=2, same binary).
# MODE 1 IS A PROBE, NOT A CANDIDATE: it keeps condensate whose heat it throws away (energy not conserved).
# PRE-REGISTERED: if the kept latent heat starves the storm track, stm_1's 35-65 column is cooler than stm_2's
# (~0.5 K) and its 35-65 precipitation is higher, toward ra_10's; if stm_1 = stm_2 in that band, the heat is NOT
# the cause and the starvation is the closure's water path (E from level 1, no re-pin). Exit 0, zero NaN.
set -u; cd "$(dirname "$0")"; rm -f STM_DONE
until [ -e RA_DONE ]; do sleep 30; done
for t in 2 1; do mkdir output_stm_$t || { echo "output_stm_$t exists"; touch STM_DONE; exit 1; }
  [ -e config_stm_$t.xml ] && { echo "config_stm_$t.xml exists"; touch STM_DONE; exit 1; }
  sed "s#output_o2val/#output_stm_$t/#" config_o2val.xml > config_stm_$t.xml; done
echo "start $(date +%H:%M)  atm_sdr $(md5sum < ../cli/atm_sdr | cut -c1-8)"
D="ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_SR_DIAG=1 ATM_PRECIP_UPWIND=1 ATM_SNOW_WINDOW=2 ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0"
run(){ ( env OMP_NUM_THREADS=6 $D $2 ../cli/atm_sdr config_stm_$1.xml > stm_$1.log 2>&1
         echo "stm_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' stm_$1.log)  $(date +%H:%M)" ) & }
run 2 ""
run 1 "ATM_RK_SCALAR_SYNC=1"
wait
for t in 2 1; do echo "== stm_$t $(grep -o 'RK_SCALAR_SYNC=[^ ]*\|WATER_CLOSURE=[^ ]*' stm_$t.log | head -2 | tr '\n' ' ')"
  grep "by |latitude|" stm_$t.log | tail -1; grep "model .*NASA .*bias" stm_$t.log | tail -1 | cut -c1-150
  grep "water budget closure" stm_$t.log | tail -1; grep 'max v-component' stm_$t.log | tail -1 | cut -c1-80
  tail -1 output_stm_$t/convergence.csv; done
touch STM_DONE
