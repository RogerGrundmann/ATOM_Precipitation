#!/bin/bash
# 2026-09-30: plan D re-run of the XML arms after making TwoCatIce::q_c_crit lazy (it was read at static init,
# before LoadConfig, and the <knobs> guard refused the config). Arms X, Y and T of run_vxml.sh; A-E there passed.
set -u; cd "$(dirname "$0")"; rm -f VXML2_DONE
. ./verify_lib.sh
./build_o0.sh WORKTREE xml2 || { touch VXML2_DONE; exit 1; }
for t in x t y; do mkdir output_vxml2_$t || { touch VXML2_DONE; exit 1; }
  sed "s#output_vxml_$t/#output_vxml2_$t/#" config_vxml_$t.xml > config_vxml2_$t.xml
  [ $t = y ] && { cp output_oc_ctl/hyd_restart_0Ma_1000.bin output_vxml2_y/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_vxml2_y/; }
done
arm vxml2_x ../cli/xml2_atm config_vxml2_x.xml
arm vxml2_y ../cli/xml2_hyd config_vxml2_y.xml
wait_arms
env OMP_NUM_THREADS=1 ../cli/xml2_atm config_vxml2_t.xml > vxml2_t.log 2>&1; te=$?
for t in x y; do [ -f output_vxml2_$t/RUN_CONFIG.txt ] && sed -i "s#output_vxml2_[a-z]/#OUT/#" output_vxml2_$t/RUN_CONFIG.txt; done
cmp_dirs vxml2_x vxml_d "X  atm new with <knobs> == old+env"
cmp_dirs vxml2_y vxml_z "D  hyd new with <knobs> == old+env"
echo "T  typo in <knobs>: exit $te; message: $(grep -a -m1 'knob' vxml2_t.log | tail -c 160)"
echo "--- banner of the XML arm:"; grep -a '\[RUN CONFIG\]' vxml2_x.log | grep -a 'RH_STORM\|+ = from' | cut -c1-220
touch VXML2_DONE
