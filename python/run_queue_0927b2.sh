#!/bin/bash
# 2026-09-27 phase B, stage 2: three lanes of ~8 threads (the model is memory-bandwidth bound: it stops gaining at
# ~8 threads, and 3 concurrent 8-thread runs each slow ~2.2x -- lanes buy ~30 % throughput, not 3x).
#   usage: bash run_queue_0927b2.sh O2   -> cli/atm_O2h / cli/hyd_O2h (4aec4dd at -O2; only if o2val PASSED)
#          bash run_queue_0927b2.sh O0   -> each script's own -O0 binary
# Every arm script still gates on its own byte check (-O0 vs -O0, valid for the default branch of 4aec4dd).
#   lane A (8 thr):  rksq_on (100) -> gpaqb (600)
#   lane B (6x1):    fx4_100 (6 arms) -> snw40 (6 arms)
#   lane C (ocean):  osf200 (3x3) -> ohs600 (2x4)
# snw600 is NOT launched here -- only if snw40 shows a large deletion (decided by reading its SNOW DIAG).
set -u; cd "$(dirname "$0")"; rm -f QUEUE_B2_DONE
case "${1:-}" in
  O2) A="ATMBIN=../cli/atm_O2h"; H="HYDBIN=../cli/hyd_O2h" ;;
  O0) A="ATMBIN_UNUSED=1";       H="HYDBIN_UNUSED=1" ;;
  *)  echo "usage: $0 O2|O0"; exit 1 ;;
esac
echo "phase B stage 2 ($1) start $(date +%H:%M)"
( env $A NT=8 bash run_rksq100.sh > run_rksq100.out 2>&1; env $A NT=8 bash run_gpaq600.sh > run_gpaq600.out 2>&1
  echo "lane A done $(date +%H:%M)" ) &
( env $A NT=1 bash run_fx4_100.sh > run_fx4_100.out 2>&1; env $A NT=1 bash run_snw40.sh > run_snw40.out 2>&1
  echo "lane B done $(date +%H:%M)" ) &
( env $H NT=3 bash run_osf200.sh > run_osf200.out 2>&1;   env $H NT=4 bash run_ohs600.sh > run_ohs600.out 2>&1
  echo "lane C done $(date +%H:%M)" ) &
wait
echo "phase B stage 2 end $(date +%H:%M)"
touch QUEUE_B2_DONE
