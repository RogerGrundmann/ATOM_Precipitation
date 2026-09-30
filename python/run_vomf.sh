#!/bin/bash
# 2026-09-30: OCN-METRIC flipped ON as defaults (user): HYD_METRIC_RADIUS 0 -> 6370, HYD_RUN_NEUMANN 0 -> 1,
# HYD_A_H_BIHARM 0 -> 1.0e19, HYD_A_H_BIHARM_SCALED 0 -> 1. Byte check both directions, -O0 both sides, 1 thread,
# restart 1000 -> 1020 from oc_ctl: old = HEAD 1b0f343 (defaults off), new = HEAD + the flip (working tree).
#   A  new clean == old + the four set;  B  new with the four = 0 == old clean;  C  new clean != old clean
#   (RUN_CONFIG may differ only by the four banner tokens)
set -u; cd "$(dirname "$0")"; rm -f VOMF_DONE
. ./verify_lib.sh
./build_o0.sh 1b0f343 omf_old  || { touch VOMF_DONE; exit 1; }
./build_o0.sh WORKTREE omf_new || { touch VOMF_DONE; exit 1; }
for t in a b c d; do mkdir output_vomf_$t || { touch VOMF_DONE; exit 1; }
  sed -e "s#output_sm_bc1/#output_vomf_$t/#" -e "s#<nm>200</nm>#<nm>20</nm>#" \
      -e "s#<checkpoint_save_iter>1200</checkpoint_save_iter>#<checkpoint_save_iter>1020</checkpoint_save_iter>#" config_sm_bc1.xml > config_vomf_$t.xml
  cp output_oc_ctl/hyd_restart_0Ma_1000.bin output_vomf_$t/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_vomf_$t/; done
ON="HYD_METRIC_RADIUS=6370 HYD_RUN_NEUMANN=1 HYD_A_H_BIHARM=1.0e19 HYD_A_H_BIHARM_SCALED=1"
OFF="HYD_METRIC_RADIUS=0 HYD_RUN_NEUMANN=0 HYD_A_H_BIHARM=0 HYD_A_H_BIHARM_SCALED=0"
arm vomf_a ../cli/omf_new_hyd config_vomf_a.xml
arm vomf_b ../cli/omf_old_hyd config_vomf_b.xml $ON
arm vomf_c ../cli/omf_new_hyd config_vomf_c.xml $OFF
arm vomf_d ../cli/omf_old_hyd config_vomf_d.xml
wait_arms
for t in a b c d; do [ -f output_vomf_$t/RUN_CONFIG.txt ] && sed -i 's#output_vomf_[a-d]/#OUT/#' output_vomf_$t/RUN_CONFIG.txt; done
cmp_dirs vomf_a vomf_b "A  FLIP == old+knobs"
cmp_dirs vomf_c vomf_d "B  new+knobs=0 == old clean"
want_differ vomf_a vomf_d "C  CONTROL new clean vs old clean -- MUST differ"
for t in a b c d; do echo "vomf_$t: NaN $(grep -c 'NaN/Inf DETECTED' vomf_$t.log 2>/dev/null) $(grep -a -o 'METRIC_RADIUS=[^ ]*  RUN_NEUMANN=[^ ]*' vomf_$t.log | head -1) $(grep -a -o 'A_H_BIHARM=[^ ]*  A_H_BIHARM_SCALED=[^ ]*' vomf_$t.log | head -1)"; done
touch VOMF_DONE
