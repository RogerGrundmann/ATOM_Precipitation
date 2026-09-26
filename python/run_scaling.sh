#!/bin/bash
# THREAD-SCALING + OPTIMISATION TEST (2026-09-26, for tomorrow; LAUNCH BY HAND on an otherwise IDLE machine --
# anything else running falsifies it).
# FOUND 2026-09-26: the Makefile has NO -O level, so every binary in this tree is built at -O0 (the binary's own
# build record: "-march=alderlake ... -ggdb -std=c++17 -ffast-math -fopenmp", no -O). Two questions:
#   (1) how much faster is -O2? cli/atm_O2 (built in a separate worktree at eb2181c with CFLAGS + -O2, Makefile
#       untouched) against cli/atm_fx4 (-O0, same source).
#   (2) with each build, does the rate stop improving past ~8 threads (memory-bandwidth / E-core bound)? That
#       decides the hardware question on the i7-13700K (8P+8E cores, 2 memory channels).
# Default configuration, 30 iterations from scratch, checkpoint every 10: the rate is taken from the momentum-budget
# file times at 10 -> 30 (20 iterations), so setup is excluded. Serial 4/8/12/16/24 threads, then 3 x 8 concurrent,
# for each build. Default thread placement, as production uses. ~1 h in total (the -O0 half ~40 min).
# (3) a first look at whether -O2 changes results: Precip / r / max|u| at iteration 30, both builds at 8 threads.
#     -O2 + -ffast-math reorders floating point, so last-digit changes are EXPECTED; a validation 600 against
#     sgzb_ctl is the real test, not this.
set -u; cd "$(dirname "$0")"; rm -f SCAL_DONE
rate(){ local d=output_$1
        local a=$(stat -c %Y $d/w_momentum_budget_10.csv 2>/dev/null) b=$(stat -c %Y $d/w_momentum_budget_30.csv 2>/dev/null)
        [ -n "$a" ] && [ -n "$b" ] && echo "scale=2; ($b - $a)/20" | bc || echo "n/a"; }
for p in scal scalo2; do for t in 4 8 12 16 24 c1 c2 c3; do
  mkdir output_${p}_$t || { echo "output_${p}_$t exists -- abort"; touch SCAL_DONE; exit 1; }; done; done
echo "start $(date +%H:%M)  load $(cut -d' ' -f1 /proc/loadavg)"
for b in "scal:atm_fx4:-O0" "scalo2:atm_O2:-O2"; do
  p=${b%%:*}; r=${b#*:}; bin=${r%%:*}; lab=${r#*:}
  for n in 4 8 12 16 24; do
    env OMP_NUM_THREADS=$n ../cli/$bin config_${p}_$n.xml > ${p}_$n.log 2>&1
    echo "$lab serial $n threads: exit $?  $(rate ${p}_$n) s/iter"
  done
  for c in c1 c2 c3; do ( env OMP_NUM_THREADS=8 ../cli/$bin config_${p}_$c.xml > ${p}_$c.log 2>&1 ) & done
  wait
  echo "$lab concurrent 3 x 8 threads: $(rate ${p}_c1) / $(rate ${p}_c2) / $(rate ${p}_c3) s/iter each"
done
echo "--- results at iteration 30, 8 threads (last-digit differences expected under -O2 -ffast-math):"
for p in scal scalo2; do
  echo "$p: $(grep 'model .*NASA' ${p}_8.log | tail -1 | sed 's/^ *//' | cut -c1-110)"
  echo "$p: $(grep 'max u-component' ${p}_8.log | tail -1 | sed 's/^ *//' | cut -c1-70)"
done
echo "end $(date +%H:%M)"
touch SCAL_DONE
