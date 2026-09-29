#!/bin/bash
# 2026-09-29: SURF-DRAG re-run. The 09-28 pair (run_sd0928.sh) did `cp -n` into an EXISTING cli/atm_sd (09-24,
# BC_SECOND_ORDER=1*), so both arms ran a stale binary. This builds 19f9a0d -O2 fresh and REFUSES an existing target.
# Pair, 600 from scratch, 2 x 8 threads, one binary, momentum budgets as in config_o2val.
#   sdr_ctl  ATM_SURF_DRAG_CONSISTENT unset (0.0)
#   sdr_1    ATM_SURF_DRAG_CONSISTENT=1.0
# PRE-REGISTERED (as 09-28): drag_sfc/coriolis below 300 m, 20-70 deg, ~1e-6 -> O(0.1-1); climate null
# (e-folding at s=1 is 4.31e5 iterations); exit 0, zero NaN; max|v| ~2.5 m/s in both (no seam mode).
set -u; cd "$(dirname "$0")"; rm -f SDR_DONE
[ -e ../cli/atm_sdr ] && { echo "cli/atm_sdr exists -- refusing"; touch SDR_DONE; exit 1; }
W=$(mktemp -d /tmp/atom_o2_XXXXXX)
(cd .. && git worktree add --detach $W 19f9a0d >/dev/null && cd $W && make -j6 atm > build.log 2>&1) \
  && cp $W/cli/atm ../cli/atm_sdr; (cd .. && git worktree remove --force $W; git worktree prune)
[ -e ../cli/atm_sdr ] || { echo "O2 build failed"; touch SDR_DONE; exit 1; }
for t in ctl 1; do mkdir output_sdr_$t || { touch SDR_DONE; exit 1; }
  [ -e config_sdr_$t.xml ] && { echo "config_sdr_$t.xml exists"; touch SDR_DONE; exit 1; }
  sed "s#output_o2val/#output_sdr_$t/#" config_o2val.xml > config_sdr_$t.xml; done
echo "start $(date +%H:%M)  atm_sdr $(md5sum < ../cli/atm_sdr | cut -c1-8)"
run(){ ( env OMP_NUM_THREADS=8 $2 ../cli/atm_sdr config_sdr_$1.xml > sdr_$1.log 2>&1
         echo "sdr_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' sdr_$1.log)  $(date +%H:%M)" ) & }
run ctl ""
run 1   "ATM_SURF_DRAG_CONSISTENT=1.0"
wait
for t in ctl 1; do echo "== sdr_$t $(grep -o 'SURF_DRAG_CONSISTENT=[^ ]*' sdr_$t.log | head -1)  $(grep -o 'BC_SECOND_ORDER=[^ ]*' sdr_$t.log | head -1)"
  echo "   max|v| every 10th print: $(grep 'max v-component' sdr_$t.log | awk 'NR%10==0{printf "%s ", $5}')"
  grep "by |latitude|" sdr_$t.log | tail -1; grep "model .*NASA .*bias" sdr_$t.log | tail -1 | cut -c1-150
  grep -E "max (v|w)-component" sdr_$t.log | tail -2 | cut -c1-100; tail -1 output_sdr_$t/convergence.csv; done
touch SDR_DONE
