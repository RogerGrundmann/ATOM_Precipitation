#!/bin/bash
# 2026-09-27 phase B, stage 1: the -O2 validation, ALONE (scaling test: -O2 is 3.4-3.9x faster than -O0).
# Waits for phase A's byte checks and the snow byte check, then run_o2val600.sh at 16 threads (~20 min).
# Stage 2 (all of phase B on cli/atm_O2h / cli/hyd_O2h if o2val passes, else -O0) is launched by hand after scoring.
set -u; cd "$(dirname "$0")"; rm -f QUEUE_B1_DONE
until [ -f QUEUE_A_DONE ] && [ -f VSNW_VERIFY_DONE ]; do sleep 30; done
echo "o2val launch $(date +%H:%M)"
NT=16 bash run_o2val600.sh > run_o2val600.out 2>&1
cat run_o2val600.out
touch QUEUE_B1_DONE
