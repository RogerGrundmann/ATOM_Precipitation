#!/bin/bash
# 2026-09-29: CONV-DEAD -- ATM_MC_ED_AREA on the energy-conserving working branch. 600 -> 640 restart from stm_2
# (closure sync 2, moisture filter + vertical T filter off, UPWIND + SNOW_WINDOW=2), fresh -O2 build cli/atm_ced,
# 2 x 8 threads, ATM_MC_DIAG on. PRE-REGISTERED: e_d share of generation falls from 99.96 % toward ~84 % or less
# (sigma^(5/9) ~ 0.16 on the demand), so P_conv at the ground rises from ~0.16 mm/a by orders of magnitude.
set -u; cd "$(dirname "$0")"; rm -f CED_DONE
[ -e ../cli/atm_ced ] && { echo "cli/atm_ced exists -- refusing"; touch CED_DONE; exit 1; }
W=$(mktemp -d /tmp/atom_o2_XXXXXX)
(cd .. && git worktree add --detach $W HEAD >/dev/null && git diff HEAD -- atmosphere > $W/.wt.diff \
   && cd $W && git apply .wt.diff && make -j6 atm > build.log 2>&1) && cp $W/cli/atm ../cli/atm_ced
(cd .. && git worktree remove --force $W; git worktree prune)
[ -e ../cli/atm_ced ] || { echo "O2 build failed"; touch CED_DONE; exit 1; }
for t in 0 1; do mkdir output_ced_$t || { touch CED_DONE; exit 1; }
  sed -e "s#output_o2val/#output_ced_$t/#" -e "s#<restart_from_iter>-1<#<restart_from_iter>600<#" -e "s#<nm>600<#<nm>640<#" \
      config_o2val.xml > config_ced_$t.xml; cp output_stm_2/atm_restart_0Ma_600.bin output_ced_$t/; done
BR="ATM_MC_DIAG=1 ATM_CWB_DIAG=1 ATM_PRECIP_UPWIND=1 ATM_SNOW_WINDOW=2 ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0"
echo "start $(date +%H:%M)  atm_ced $(md5sum < ../cli/atm_ced | cut -c1-8)"
for t in 0 1; do ( env OMP_NUM_THREADS=8 $BR ATM_MC_ED_AREA=$t ../cli/atm_ced config_ced_$t.xml > ced_$t.log 2>&1
  echo "ced_$t exit $?  NaN $(grep -c 'NaN/Inf DETECTED' ced_$t.log)  $(date +%H:%M)" ) & done
wait
for t in 0 1; do echo "== ced_$t"; grep -a '\[MC DIAG\] column budget\|\[MC DIAG\] as shares\|\[MC DIAG\] active' ced_$t.log | tail -3 | cut -c1-200
  grep -a 'P_conv mean' ced_$t.log | tail -1; grep -a "model .*NASA .*bias" ced_$t.log | tail -1 | cut -c1-120
  grep -a "by |latitude|" ced_$t.log | tail -1; grep -a 'max v-component\|MC_t.*truncated\|max temperature' ced_$t.log | tail -2 | cut -c1-90; done
touch CED_DONE
