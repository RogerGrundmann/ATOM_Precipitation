#!/bin/bash
# 2026-09-29: RAIN-CONV -- re-sweep ATM_RAIN_AREA on the CONVERGED rain branch. 0.10 was fitted 2026-09-01 on the
# shipped 3-pass rain column, which never converges (the NASA match was pass 3 of an oscillation). On the upwind
# branch the column converges and rains 85.9 mm/a (rp_uw), 173.3 with ATM_SNOW_WINDOW=2 (rp_uw2); S_ev removes ~91 %
# of the rain sources. S_ev's grid-mean tendency goes as f^(5/9), so a smaller shaft fraction evaporates less.
# Branch: ATM_PRECIP_UPWIND=1 ATM_SNOW_WINDOW=2. Binary cli/atm_sdr (fresh -O2 of 19f9a0d, md5 7fee2012).
# 4 arms, 600 from scratch, 4 x 6 threads, concurrent; starts after run_sdr.sh (SDR_DONE).
#   ra_10   RAIN_AREA 0.10  (control; must reproduce rp_uw2's 173 mm/a to the thread noise)
#   ra_03   0.03
#   ra_01   0.01
#   ra_003  0.003
# PRE-REGISTERED: (1) precipitation rises monotonically as f falls (the 09-01 non-monotonic peak came from the floor
# injection, which ATM_ICE_LIMIT_ARRIVING removed); (2) with S_ev at 91 % of sources, even S_ev -> 0 gives at most
# ~10x, i.e. ~1000-1700 mm/a -- if 0.003 is still far below NASA, the gap is not evaporation but generation;
# (3) judged on r, sigma, the four bands and land/ocean, NOT the mean; (4) exit 0, zero NaN.
# NOTE: f below ~0.05 is not a physical stratiform shaft fraction; any value picked here is a scaffold constant.
set -u; cd "$(dirname "$0")"; rm -f RA_DONE
until [ -e SDR_DONE ]; do sleep 30; done
T=(10 03 01 003); F=(0.10 0.03 0.01 0.003)
for t in "${T[@]}"; do mkdir output_ra_$t || { echo "output_ra_$t exists"; touch RA_DONE; exit 1; }
  [ -e config_ra_$t.xml ] && { echo "config_ra_$t.xml exists"; touch RA_DONE; exit 1; }
  sed "s#output_o2val/#output_ra_$t/#" config_o2val.xml > config_ra_$t.xml; done
echo "start $(date +%H:%M)  atm_sdr $(md5sum < ../cli/atm_sdr | cut -c1-8)"
D="ATM_CWB_DIAG=1 ATM_SR_DIAG=1 ATM_RAIN_PASS_DIAG=1 ATM_PRECIP_UPWIND=1 ATM_SNOW_WINDOW=2"
for n in 0 1 2 3; do t=${T[$n]}
  ( env OMP_NUM_THREADS=6 $D ATM_RAIN_AREA=${F[$n]} ../cli/atm_sdr config_ra_$t.xml > ra_$t.log 2>&1
    echo "ra_$t exit $?  NaN $(grep -c 'NaN/Inf DETECTED' ra_$t.log)  $(date +%H:%M)" ) & done
wait
for t in "${T[@]}"; do echo "== ra_$t $(grep -o 'RAIN_AREA=[^ ]*' ra_$t.log | head -1)"
  grep "by |latitude|" ra_$t.log | tail -1; grep "model .*NASA .*bias" ra_$t.log | tail -1 | cut -c1-150
  grep -i "land .*ocean" ra_$t.log | tail -1 | cut -c1-120; grep "SR DIAG\] as shares" ra_$t.log | tail -1
  grep 'max v-component' ra_$t.log | tail -1 | cut -c1-80; tail -1 output_ra_$t/convergence.csv; done
touch RA_DONE
