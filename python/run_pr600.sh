#!/bin/bash
# ATM_RH_MIN_PTOP RE-SWEEP ON THE 2026-09-24 DEFAULT: ATM_RAD_TOPO=1 (flip) + second-order BCs + PTOP 490.
# WHY: rt600 (RAD_TOPO=1, PTOP 475) raised the parity-mean precipitation +5.5 % over el2 (907.1 -> 957.0),
# so PTOP 490, fitted without RAD_TOPO (965.2 = -1.3 %), should now overshoot NASA.
# Binary cli/atm_rt (post-flip, md5 4305c587). 600 from scratch, 4 arms x 6 threads concurrent.
# Waits for the hsf 600->1200 pair (HSF1200_DONE) AND the RAD_TOPO byte check (RADTOPO_VERIFY_DONE).
# pr490 is ALSO the new default's own 600-iteration arm (compare hsf_ctl: same, RAD_TOPO=0).
# PRE-REGISTERED (from the 2026-09-23 sweep scaled by rt600's +5.5 %):
#   1. parity-mean precip ~ 915 / 957 / 985 / 1018 for 465 / 475 / 482 / 490; NASA 978.3 lands at ~480.
#   2. PTOP stays a pure scale knob: 35-65 ~205, 65-90 ~24 flat to ~1 %; sigma/mean flat.
#   3. land/ocean RATIO stays ~0.77 (the RAD_TOPO shape change) at every PTOP.
#   4. r ~0.47-0.48, no real optimum; no drift 200 -> 600.
#   5. clean: exit 0, zero NaN through 155/357/483; pr490 max|v| seam mode as in hsf_ctl (~26 m/s).
set -u; cd "$(dirname "$0")"; rm -f PR600_DONE
until [ -f HSF1200_DONE ] && [ -f RADTOPO_VERIFY_DONE ]; do sleep 60; done
echo "start $(date +%H:%M)"
for V in 465 475 482 490; do
  mkdir -p output_pr$V
  ( env OMP_NUM_THREADS=6 ATM_MC_CAP_DIAG=1 ATM_RH_MIN_PTOP=$V \
        ../cli/atm_rt config_pr$V.xml > pr$V.log 2>&1
    echo "pr$V exit $?  NaN $(grep -c 'NaN/Inf DETECTED' pr$V.log)  $(date +%H:%M)" ) &
done
wait; touch PR600_DONE
