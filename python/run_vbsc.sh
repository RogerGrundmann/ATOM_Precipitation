#!/bin/bash
# 2026-09-29: HYD_A_H_BIHARM_SCALED (new, default 0) -- -O0 off-branch byte check, 1 thread, restart 1000 -> 1020
# from oc_ctl. old = cli/sqr_new_hyd (113221f -O0), new = cli/bsc_hyd (113221f + the knob, -O0).
#   A  new clean == old clean;  B  new B=3e17 == old B=3e17 (knob unset, biharmonic ON);
#   C  new B=3e17 SCALED=1 != old B=3e17  (control: the knob must fire)
# B = 3e17 is below the ~9e17 explicit limit at floor 0.26, so every arm is stable.
set -u; cd "$(dirname "$0")"; rm -f VBSC_DONE
. ./verify_lib.sh
for t in a b c d e; do mkdir output_vbsc_$t || { touch VBSC_DONE; exit 1; }
  sed -e "s#output_sm_bc1/#output_vbsc_$t/#" -e "s#<nm>200</nm>#<nm>20</nm>#" \
      -e "s#<checkpoint_save_iter>1200</checkpoint_save_iter>#<checkpoint_save_iter>1020</checkpoint_save_iter>#" config_sm_bc1.xml > config_vbsc_$t.xml
  cp output_oc_ctl/hyd_restart_0Ma_1000.bin output_vbsc_$t/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_vbsc_$t/; done
arm vbsc_a ../cli/bsc_hyd     config_vbsc_a.xml
arm vbsc_b ../cli/sqr_new_hyd config_vbsc_b.xml
arm vbsc_c ../cli/bsc_hyd     config_vbsc_c.xml HYD_A_H_BIHARM=3e17
arm vbsc_d ../cli/sqr_new_hyd config_vbsc_d.xml HYD_A_H_BIHARM=3e17
arm vbsc_e ../cli/bsc_hyd     config_vbsc_e.xml HYD_A_H_BIHARM=3e17 HYD_A_H_BIHARM_SCALED=1
wait_arms
for t in a b c d e; do [ -f output_vbsc_$t/RUN_CONFIG.txt ] && sed -i 's#output_vbsc_[a-e]/#OUT/#' output_vbsc_$t/RUN_CONFIG.txt; done
cmp_dirs vbsc_a vbsc_b "A  new clean == old clean"
cmp_dirs vbsc_c vbsc_d "B  new B=3e17 == old B=3e17"
want_differ vbsc_e vbsc_d "C  CONTROL new B=3e17 SCALED=1 vs old B=3e17 -- MUST differ"
touch VBSC_DONE
