#!/bin/bash
# 2026-09-28: SURF-DRAG (B.9) -- the surface Rayleigh drag is 4.0e5 too weak (a spurious *dt and L_atm instead
# of metricShellLength()). ATM_SURF_DRAG_CONSISTENT=<s> blends to the consistent coefficient; never run.
# Pair, 600 from scratch, today's defaults, -O2 build of 9438b0c (cli/atm_sd), 2 x 8 threads, momentum budgets on
# (written every checkpoint). Starts after run_om0928.sh (user: "after 15:15").
#   sd_ctl  ATM_SURF_DRAG_CONSISTENT unset (0.0)
#   sd_1    ATM_SURF_DRAG_CONSISTENT=1.0
# PRE-REGISTERED: (1) connection: w_momentum_budget drag_sfc/coriolis in the lowest 300 m rises from ~1e-6 to
# O(0.1) (the coefficient ratio drag/force_nd = kf/omega = 0.159). (2) The e-folding at s=1 is 4.31e5 iterations,
# so 600 iterations are 0.14 % of one: surface winds move <~1 %, climate (precip, r, bands) null. (3) exit 0,
# zero NaN. A large change would mean the drag acts through something faster than its own e-folding.
set -u; cd "$(dirname "$0")"; rm -f SD0928_DONE
until [ -e OM0928_DONE ] && [ "$(date +%H%M)" -ge 1515 ]; do sleep 120; done
W=$(mktemp -d /tmp/atom_o2_XXXXXX)
(cd .. && git worktree add --detach $W 9438b0c >/dev/null && cd $W && make -j6 atm > build.log 2>&1) \
  && cp -n $W/cli/atm ../cli/atm_sd; (cd .. && git worktree remove --force $W; git worktree prune)
[ -e ../cli/atm_sd ] || { echo "O2 build failed"; touch SD0928_DONE; exit 1; }
for t in ctl 1; do mkdir output_sd_$t || { touch SD0928_DONE; exit 1; }
  sed "s#output_o2val/#output_sd_$t/#" config_o2val.xml > config_sd_$t.xml; done
echo "start $(date +%H:%M)  atm_sd $(md5sum < ../cli/atm_sd | cut -c1-8)"
run(){ ( env OMP_NUM_THREADS=${NT:-8} $2 ../cli/atm_sd config_sd_$1.xml > sd_$1.log 2>&1
         echo "sd_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' sd_$1.log)  $(date +%H:%M)" ) & }
run ctl ""
run 1   "ATM_SURF_DRAG_CONSISTENT=1.0"
wait
for t in ctl 1; do echo "== sd_$t $(grep -o 'SURF_DRAG_CONSISTENT=[^ ]*' sd_$t.log | head -1)"; grep "by |latitude|" sd_$t.log | tail -1
  grep "model .*NASA .*bias" sd_$t.log | tail -1 | cut -c1-150; grep -E "max (v|w)-component" sd_$t.log | tail -2 | cut -c1-100
  tail -1 output_sd_$t/convergence.csv; done
touch SD0928_DONE
