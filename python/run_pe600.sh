#!/bin/bash
# ATM_RH_MIN_PTOP RE-SWEEP ON THE 2026-09-22 DEFAULT (ATM_MC_EVAP_LIMIT=1 + ATM_MC_T_NDIM=1.0).
# WHY: PTOP=475 was fitted on the branch whose precipitation drifted upward through unbacked
# sub-cloud evaporation moistening via MC_q. With the limiter on, the drift is gone and the
# mean at 600 is 927.8 mm/a (high parity), -4.5 % of NASA. Last re-sweep (2026-09-21) found
# the bands INVARIANT to PTOP and sigma tracking the mean to 0.3 %, so the discriminators are
# r and the land/ocean partition (el2: 706.5 / 1015.4 against NASA 782.3 / 1055.8).
# FREE ARM: output_el2 IS PTOP=475 on this binary (cli/atm_el, md5 d5116dd7), same env.
# DIRECTION: UPWARD this time (mean is below NASA); 450 is the direction control.
# PRE-REGISTERED:
#   1. mean ~ +100 mm/a per 25 hPa, monotone (2026-09-21: +104.7 for 475->500). So 490 ~ 990,
#      500 ~ 1030, 525 ~ 1130, 450 ~ 830.
#   2. 35-65 / 65-90 bands flat to ~1 % across all arms (184.1 / 20.5 at 475).
#   3. sigma rises with the mean in proportion (pure scaling); if sigma/mean is NOT flat the
#      limiter changed how the floor acts and the old sweep does not transfer.
#   4. land/ocean: both rise; land needs +11 %, ocean +4 %, so NO single PTOP fixes both.
#   5. drift: el2's precip was flat 100->600; each arm's trajectory must stay flat too, or PTOP
#      re-opens the drift.
# THREADS: 4 arms x 6 = 24 concurrent. Measure s/iter from checkpoint mtimes.
set -u; cd "$(dirname "$0")"; rm -f PE600_DONE
for V in 450 490 500 525; do
  mkdir -p output_pe$V
  ( env OMP_NUM_THREADS=6 ATM_MC_CAP_DIAG=1 ATM_MC_EVAP_LIMIT=1 ATM_MC_T_NDIM=1 ATM_RH_MIN_PTOP=$V \
        ../cli/atm_el config_pe$V.xml > pe$V.log 2>&1
    echo "pe$V exit $?  NaN $(grep -c 'NaN/Inf DETECTED' pe$V.log)  $(date +%H:%M)" ) &
done
wait; touch PE600_DONE
