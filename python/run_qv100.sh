#!/bin/bash
# TEST of the closure-runaway mechanism (2026-09-25): is the moisture filter's VERTICAL pass the engine?
# 100 from scratch, cli/atm_dqv, 2 arms x 3 threads, alongside sgzb600. CWB + MC_CAP diag on, as bis_*.
#   qv_on   ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0   -- compare bis_ctl (closure on, filter vertical ON)
#   qv_off  ATM_DAMP_Q_VERT=0                       -- compare bis_wc  (closure off): the filter's effect
#                                                      on the shipped branch, where RK4 discards it
# PRE-REGISTERED: qv_on's climb over 40-100 falls from bis_ctl's 11.87 mm/a/iter to ~bis_wc's 3.6 or less;
# 15-35 / 65-90 bands at 100 near bis_wc (217 / 16), not bis_ctl (1488 / 350); CWB damp_wiggles(q) surface/
# aloft transfer ~0; the evaporation row turns POSITIVE (skin copy no longer mirrors a drained level 1);
# microphysics sink near the closure-off level. qv_off: ~null against bis_wc (the pump is discarded there).
set -u; cd "$(dirname "$0")"; rm -f QV100_DONE
for t in on off; do mkdir output_qv_$t || exit 1; done
echo "start $(date +%H:%M)"
run(){ local tag=$1; shift
       ( env OMP_NUM_THREADS=3 ATM_CWB_DIAG=1 ATM_MC_CAP_DIAG=1 "$@" ../cli/atm_dqv config_$tag.xml > $tag.log 2>&1
         echo "$tag exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $tag.log)  $(date +%H:%M)" ) & }
run qv_on  ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0
run qv_off ATM_DAMP_Q_VERT=0
wait
touch QV100_DONE
