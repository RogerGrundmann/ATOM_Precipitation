#!/bin/bash
# 2026-09-30: ATM_HYDRO_SPLIT 1200 -> 1800 continuation PAIR from each arm's own iteration-1200 checkpoint (user: keep the
# old checkpoints). SAME binary cli/atm_hsf (-O0, 09-23, built from an uncommitted tree, so it cannot be rebuilt at -O2
# faithfully) and the same env as run_hsf1200.sh -- the only change is the length. Queued behind rhst_115; 6 threads
# each because the om_b* ocean pair holds 12. nm is the TOTAL iteration count in the atmosphere.
# At 1200 (hsf1200_on): max|p_dyn| pre-clip 0.192 (600) -> 0.265 (1200), growth slowing ~35 %; ceiling 3.0.
# PRE-REGISTERED:
#   1. max|p_dyn| growth per 100 iterations keeps falling (plateau) and stays far below 3.0, zero CEILING BINDING
#      -> ATM-HSPLIT closes. Stays linear or re-accelerates -> investigate the projection.
#   2. max|u| tracks ctl; max|v| seam mode not worse than ctl. 3. climate null vs ctl at 1800. 4. exit 0, zero NaN.
set -u; cd "$(dirname "$0")"; rm -f HSF1800_DONE
until [ -e RHST115_DONE ]; do sleep 60; done
for a in ctl on; do
  mkdir output_hsf1800_$a || { touch HSF1800_DONE; exit 1; }
  [ -e config_hsf1800_$a.xml ] && { touch HSF1800_DONE; exit 1; }
  sed -e "s#output_hsf1200_$a/#output_hsf1800_$a/#" -e 's#<checkpoint_save_iter>1200<#<checkpoint_save_iter>1800<#' \
      -e 's#<restart_from_iter>600<#<restart_from_iter>1200<#' -e 's#<nm>1200<#<nm>1800<#' \
      config_hsf1200_$a.xml > config_hsf1800_$a.xml
  cp -n output_hsf1200_$a/atm_restart_0Ma_1200.bin output_hsf1800_$a/
done
echo "start $(date +%H:%M)  atm_hsf $(md5sum < ../cli/atm_hsf | cut -c1-8)"
( env OMP_NUM_THREADS=6 ATM_RAD_TOPO=0 ATM_UBUD_BALANCE=1 \
      ../cli/atm_hsf config_hsf1800_ctl.xml > hsf1800_ctl.log 2>&1; echo "ctl exit $?" ) &
( env OMP_NUM_THREADS=6 ATM_RAD_TOPO=0 ATM_UBUD_BALANCE=1 ATM_HYDRO_SPLIT=1.0 \
      ../cli/atm_hsf config_hsf1800_on.xml  > hsf1800_on.log  2>&1; echo "on exit $?" ) &
wait
for a in ctl on; do
  echo "== $a  NaN $(grep -c 'NaN/Inf DETECTED' hsf1800_$a.log)  BINDING $(grep -c 'CEILING BINDING' hsf1800_$a.log)"
  echo "   max|p_dyn| pre-clip every 100: $(grep -a 'max|p_dyn| pre-clip' hsf1800_$a.log | awk 'NR%25==0{printf "%s ", $10}')"
  grep -a 'max u-component' hsf1800_$a.log | tail -1 | cut -c1-80; grep -a 'max v-component' hsf1800_$a.log | tail -1 | cut -c1-80
  grep -a "by |latitude|" hsf1800_$a.log | grep -v MFC | tail -1; tail -1 output_hsf1800_$a/convergence.csv
done
echo "$(date +%H:%M)"; touch HSF1800_DONE
# STOPPED 2026-09-30 ~09:33 at the user's instruction, at ~iteration 1330 (max|p_dyn| 0.282 in _on); partial output kept.
