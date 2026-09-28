#!/bin/bash
# 2026-09-28: HYD_METRIC_SIN_FLOOR default 0.26 (9438b0c) -- -O0 byte check both directions, 1 thread,
# restart 1000 -> 1020 from oc_ctl (shipped metric). old = cli/hsp1_new_hyd (d5d0272 -O0, floor 0.4), new = 9438b0c.
#   A  new clean == old FLOOR=0.26;  B  new FLOOR=0.4 == old clean;  C  new clean != old clean
set -u; cd "$(dirname "$0")"; rm -f VOSF0928_DONE
. ./verify_lib.sh
[ -e ../cli/osf_new_hyd ] || ./build_o0.sh 9438b0c osf_new || { touch VOSF0928_DONE; exit 1; }
for t in a b c d; do mkdir output_vosf_$t || { touch VOSF0928_DONE; exit 1; }
  sed -e "s#output_sm_bc1/#output_vosf_$t/#" -e "s#<nm>200</nm>#<nm>20</nm>#" \
      -e "s#<checkpoint_save_iter>1200</checkpoint_save_iter>#<checkpoint_save_iter>1020</checkpoint_save_iter>#" config_sm_bc1.xml > config_vosf_$t.xml
  cp output_oc_ctl/hyd_restart_0Ma_1000.bin output_vosf_$t/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_vosf_$t/; done
arm vosf_a ../cli/osf_new_hyd  config_vosf_a.xml
arm vosf_b ../cli/hsp1_new_hyd config_vosf_b.xml HYD_METRIC_SIN_FLOOR=0.26
arm vosf_c ../cli/osf_new_hyd  config_vosf_c.xml HYD_METRIC_SIN_FLOOR=0.4
arm vosf_d ../cli/hsp1_new_hyd config_vosf_d.xml
wait_arms
for t in a b c d; do [ -f output_vosf_$t/RUN_CONFIG.txt ] && sed -i 's#output_vosf_[a-d]/#OUT/#' output_vosf_$t/RUN_CONFIG.txt; done
cmp_dirs vosf_a vosf_b "A  FLIP == old+knob"
cmp_dirs vosf_c vosf_d "B  new+knob=0.4 == old clean"
want_differ vosf_a vosf_d "C  CONTROL new clean vs old clean -- MUST differ"
touch VOSF0928_DONE
