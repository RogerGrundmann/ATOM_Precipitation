#!/bin/bash
# ATM_SATADJ_FREEZE_LATENT, the full-strength RESTART pair. Chained behind SF600_DONE.
#
# Seed: output_mn_ctl/atm_restart_0Ma_600.bin (md5 082f884e) -- 600 from scratch on the CURRENT
# post-revert defaults, written by B+1's own control, so it is the cleanest seed in the tree.
# 600 -> 700 is 100 iterations: past the ~40 below which a restart prints ubud_* as all zeros.
#
# A RESTART rather than from scratch because this knob acts on cold cells that HOLD LIQUID, and
# on the slice only 0.05 % of the cloud sits below t_00 -- that population exists in a spun-up
# field and is thin in the analytic initial one. The byte check's want_differ control already
# proved the knob fires from scratch at nm = 20 (6 of 14 files differ, including BOTH momentum
# budgets), so what this pair adds is magnitude, not connectivity.
#
# ★ BUILT-IN CONTROL: ATM_CWB_DIAG's SaturationAdjust row and ATM_SATADJ_DIAG's buckets MUST NOT
# MOVE. This knob is a pure ENERGY repair -- it changes t_row[k] and no water at all. If the water
# budget shifts, the change is touching something it should not, and that is a failure, not a
# result. ATM_SATADJ_DIAG's cold_delete bucket in particular must stay exactly 0.
#
# PRE-REGISTERED: 0.3085 K per g/kg frozen, and the cold cloud is thin, so expect a LOCAL warming
# of order 1e-3 K and a null on every global metric. The interesting direction is the feedback:
# a warmer cell has a larger q_sat, so it holds more vapour and makes less ice -- i.e. this
# CORRECT repair should REDUCE cirrus slightly, against ATM_ICE_COLD which was written to increase
# it. If IWP rises instead, the sign reasoning here is wrong and that is worth recording.
set -u
cd "$(dirname "$0")"
export OMP_NUM_THREADS=12
while [ ! -f SF600_DONE ]; do sleep 60; done
rm -f FL700_DONE
echo "=== fl700 600->700 restart pair, 12 threads each $(date +%H:%M:%S)"
( env ATM_CWB_DIAG=1 ATM_SATADJ_DIAG=1                            ../cli/atm config_fl_ctl.xml > fl_ctl.log 2>&1
  echo "  fl_ctl exit $?  nan $(grep -ci 'nan' fl_ctl.log)  $(date +%H:%M:%S)" ) &
( env ATM_CWB_DIAG=1 ATM_SATADJ_DIAG=1 ATM_SATADJ_FREEZE_LATENT=1 ../cli/atm config_fl_on.xml  > fl_on.log  2>&1
  echo "  fl_on  exit $?  nan $(grep -ci 'nan' fl_on.log)   $(date +%H:%M:%S)" ) &
wait
echo "=== fl700 done $(date +%H:%M:%S)"
touch FL700_DONE
