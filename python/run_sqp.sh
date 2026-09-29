#!/bin/bash
# 2026-09-29: WHICH 09-28 FLIP BROUGHT THE SEAM MODE BACK? sd_ctl (9438b0c, SEAM_Q_CONSERVE=2 default, 600 from
# scratch) ended with max|v| 26.02 m/s at 11N 1E 236 m -- the ATM_BC_SECOND_ORDER signature -- while rp_ctl/rp_oqm
# (c1aae23 + OROG) were clean at 2.3-2.5. Suspects: ATM_SEAM_Q_CONSERVE=2 (only 20-iter tested), else
# ATM_TURB_SIN_FLOOR=1 (100-iter tested). Pair, 600 from scratch, -O2 build of 113221f (default SEAM_Q=0),
# 2 x 8 threads, one binary. Starts only if run_vsqr2.sh (113221f's byte check) PASSED.
#   sqp_0  default (ATM_SEAM_Q_CONSERVE=0)
#   sqp_2  ATM_SEAM_Q_CONSERVE=2
# PRE-REGISTERED: if SEAM_Q=2 is the cause, sqp_2 reproduces sd_ctl (max|v| grows ~1.28x per checkpoint from
# ~iter 300 to ~26 at 11N 1E) and sqp_0 stays at 2.3-2.5 m/s at 18S 68W 12 km. If BOTH grow, the cause is not
# SEAM_Q (next suspect ATM_TURB_SIN_FLOOR). Climate scores are not the test: the mode is invisible to them.
set -u; cd "$(dirname "$0")"; rm -f SQP_DONE
until [ -e VSQR2_DONE ]; do sleep 30; done
if [ "$(grep 'DIFFERS:' run_vsqr2.out | grep -vc 'RUN_CONFIG.txt')" != 0 ] || ! grep -q 'C  CONTROL.*PASS' run_vsqr2.out \
   || grep -q 'exit [^0]' run_vsqr2.out; then echo "byte check did not pass -- NOT started"; touch SQP_DONE; exit 1; fi
[ -e ../cli/atm_sqp ] && { echo "cli/atm_sqp exists"; touch SQP_DONE; exit 1; }
W=$(mktemp -d /tmp/atom_o2_XXXXXX)
(cd .. && git worktree add --detach $W 113221f >/dev/null && cd $W && make -j6 atm > build.log 2>&1) \
  && cp -n $W/cli/atm ../cli/atm_sqp; (cd .. && git worktree remove --force $W; git worktree prune)
[ -e ../cli/atm_sqp ] || { echo "O2 build failed"; touch SQP_DONE; exit 1; }
for t in 0 2; do mkdir output_sqp_$t || { touch SQP_DONE; exit 1; }
  sed "s#output_o2val/#output_sqp_$t/#" config_o2val.xml > config_sqp_$t.xml; done
echo "start $(date +%H:%M)  atm_sqp $(md5sum < ../cli/atm_sqp | cut -c1-8)"
run(){ ( env OMP_NUM_THREADS=8 $2 ../cli/atm_sqp config_sqp_$1.xml > sqp_$1.log 2>&1
         echo "sqp_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' sqp_$1.log)  $(date +%H:%M)" ) & }
run 0 ""
run 2 "ATM_SEAM_Q_CONSERVE=2"
wait
for t in 0 2; do echo "== sqp_$t $(grep -o 'SEAM_Q_CONSERVE=[^ ]*' sqp_$t.log | head -1)"
  echo "   max|v| trajectory: $(grep 'max v-component' sqp_$t.log | awk '{printf "%.2f ", $4}')"
  grep -E "max (v|w)-component" sqp_$t.log | tail -2 | cut -c1-100
  grep "by |latitude|" sqp_$t.log | tail -1; grep "model .*NASA .*bias" sqp_$t.log | tail -1 | cut -c1-150
  tail -1 output_sqp_$t/convergence.csv; done
touch SQP_DONE
