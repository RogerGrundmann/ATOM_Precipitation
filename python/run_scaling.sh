#!/bin/bash
# THREAD-SCALING TEST (2026-09-26, for tomorrow; LAUNCH BY HAND on an otherwise IDLE machine -- anything else
# running falsifies it). Question: is this model memory-bandwidth-bound on the i7-13700K (8P+8E cores, 24 threads,
# 2 memory channels)? That decides whether a faster desktop can give 2x (only with more memory channels).
# Default configuration, cli/atm_fx4, 30 iterations from scratch, checkpoint every 10: the rate is taken from the
# momentum-budget file times at 10 -> 30 (20 iterations), so setup is excluded.
#   serial: 4, 8, 12, 16, 24 threads, one at a time
#   concurrent: 3 runs x 8 threads at once (the throughput mode used for independent arms)
# READING: if s/iter stops improving past ~8 threads, the limit is memory bandwidth (and E-cores); throughput
# per day is then best with several runs in parallel, and a 2x machine needs 4-8 memory channels.
# Default thread placement, as production runs use it (no OMP_PROC_BIND), so the numbers apply to real runs.
set -u; cd "$(dirname "$0")"; rm -f SCAL_DONE
rate(){ local d=output_scal_$1
        local a=$(stat -c %Y $d/w_momentum_budget_10.csv 2>/dev/null) b=$(stat -c %Y $d/w_momentum_budget_30.csv 2>/dev/null)
        [ -n "$a" ] && [ -n "$b" ] && echo "scale=2; ($b - $a)/20" | bc || echo "n/a"; }
for t in 4 8 12 16 24 c1 c2 c3; do mkdir output_scal_$t || { echo "output_scal_$t exists -- abort"; touch SCAL_DONE; exit 1; }; done
echo "start $(date +%H:%M)  load $(cut -d' ' -f1 /proc/loadavg)"
for n in 4 8 12 16 24; do
  env OMP_NUM_THREADS=$n ../cli/atm_fx4 config_scal_$n.xml > scal_$n.log 2>&1
  echo "serial  $n threads: exit $?  $(rate $n) s/iter"
done
for c in c1 c2 c3; do ( env OMP_NUM_THREADS=8 ../cli/atm_fx4 config_scal_$c.xml > scal_$c.log 2>&1 ) & done
wait
echo "concurrent 3 x 8 threads: $(rate c1) / $(rate c2) / $(rate c3) s/iter each"
echo "end $(date +%H:%M)"
touch SCAL_DONE
