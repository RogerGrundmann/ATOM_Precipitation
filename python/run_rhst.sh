#!/bin/bash
# 2026-09-29: STORM fill hypothesis. Working branch (python/working_branch.env, now incl. ATM_SEAM_Q_CONSERVE=2),
# 600 from scratch, fresh -O2 build (cli/atm_rhst = 2182bc4 + ATM_RH_STORM), 2 x 8 threads, CWB bands + MFC on.
#   rhst_ctl    working branch
#   rhst_moist  + ATM_RH_STORM=1.25 (initial RH x1.25 at 55 deg, Gaussian width 15 deg, cap 0.98; ocean surface 0.75 -> 0.94)
# PRE-REGISTERED: if the storm-track deficit is spin-up (the column filling at +2471 mm/a), rhst_moist rains
# markedly more in 35-65 (and 65-90) from early on, and its 35-65 column gain (CWB-BANDS NET) is smaller; tropics
# unchanged. If 35-65 rain barely moves, the band is limited by condensation/re-evaporation, not by its water content.
# rhst_ctl vs cvs_all measures the SEAM_Q_CONSERVE=2 addition (35-65 seam bucket -155 mm/a -> ~0 expected).
set -u; cd "$(dirname "$0")"; rm -f RHST_DONE
[ -e ../cli/atm_rhst ] && { echo "cli/atm_rhst exists -- refusing"; touch RHST_DONE; exit 1; }
W=$(mktemp -d /tmp/atom_o2_XXXXXX)
(cd .. && git worktree add --detach $W HEAD >/dev/null && git diff HEAD -- atmosphere > $W/.wt.diff \
   && cd $W && git apply .wt.diff && make -j6 atm > build.log 2>&1) && cp $W/cli/atm ../cli/atm_rhst
(cd .. && git worktree remove --force $W; git worktree prune)
[ -e ../cli/atm_rhst ] || { echo "O2 build failed"; touch RHST_DONE; exit 1; }
for t in ctl moist; do mkdir output_rhst_$t || { touch RHST_DONE; exit 1; }
  [ -e config_rhst_$t.xml ] && { touch RHST_DONE; exit 1; }
  sed "s#output_o2val/#output_rhst_$t/#" config_o2val.xml > config_rhst_$t.xml; done
. ./working_branch.env
D="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_MFC_DIAG=1 ATM_SR_DIAG=1"
echo "start $(date +%H:%M)  atm_rhst $(md5sum < ../cli/atm_rhst | cut -c1-8)"
run(){ ( env OMP_NUM_THREADS=8 $D $2 ../cli/atm_rhst config_rhst_$1.xml > rhst_$1.log 2>&1
         echo "rhst_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' rhst_$1.log)  $(date +%H:%M)" ) & }
run ctl   ""
run moist "ATM_RH_STORM=1.25"
wait
for t in ctl moist; do echo "== rhst_$t"; grep -a "by |latitude|" rhst_$t.log | grep -v MFC | tail -1; grep -a "model .*NASA .*bias" rhst_$t.log | tail -1 | cut -c1-150
  grep -a -i "land .*ocean" rhst_$t.log | tail -1 | cut -c1-100; grep -a 'CWB-BANDS\] \(NET\|BC:phi\|evaporation\)' rhst_$t.log | tail -3
  grep -a 'water budget closure' rhst_$t.log | tail -1; grep -a 'max v-component' rhst_$t.log | tail -1 | cut -c1-80; tail -1 output_rhst_$t/convergence.csv; done
touch RHST_DONE
