#!/bin/bash
# 2026-09-30: KNOB-INV plan B -- lib/Knobs.h registry; every getenv site converted to knob:: accessors, banners
# generated. NOTHING may change. -O0 both sides, 1 thread. old = c90630c, new = c90630c + the registry (working tree).
#   atm (nm 20 from scratch):  A new clean == old clean;  B new + working_branch.env == old + working_branch.env
#   hyd (1000 -> 1020 from oc_ctl): C new clean == old clean;  D new + 5 knobs set == old + the same
#   controls: E atm new+env vs old clean MUST differ;  F hyd new+knobs vs old clean MUST differ
#   RUN_CONFIG.txt differs by construction (new banner); it is the only file allowed to.
set -u; cd "$(dirname "$0")"; rm -f VKNOB_DONE
. ./verify_lib.sh
[ -e ../cli/knob_old_hyd ] || ./build_o0.sh c90630c knob_old || { touch VKNOB_DONE; exit 1; }
./build_o0.sh WORKTREE knob_new || { touch VKNOB_DONE; exit 1; }
WB=$(grep '^export ' working_branch.env | sed 's/^export //; s/ *#.*$//' | tr '\n' ' ')
HK="HYD_METRIC_RADIUS=0 HYD_RUN_NEUMANN=0 HYD_A_H_BIHARM=0 HYD_A_H_BIHARM_SCALED=0 HYD_T_FREEZE=0"
echo "working branch env: $WB"
for t in a b c d; do mkdir output_vknob_$t || { touch VKNOB_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vknob_$t/#" config_vconv_new.xml > config_vknob_$t.xml; done
for t in e f g h; do mkdir output_vknob_$t || { touch VKNOB_DONE; exit 1; }
  sed -e "s#output_sm_bc1/#output_vknob_$t/#" -e "s#<nm>200</nm>#<nm>20</nm>#" \
      -e "s#<checkpoint_save_iter>1200</checkpoint_save_iter>#<checkpoint_save_iter>1020</checkpoint_save_iter>#" config_sm_bc1.xml > config_vknob_$t.xml
  cp output_oc_ctl/hyd_restart_0Ma_1000.bin output_vknob_$t/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_vknob_$t/; done
arm vknob_a ../cli/knob_new_atm config_vknob_a.xml
arm vknob_b ../cli/knob_old_atm config_vknob_b.xml
arm vknob_c ../cli/knob_new_atm config_vknob_c.xml $WB
arm vknob_d ../cli/knob_old_atm config_vknob_d.xml $WB
arm vknob_e ../cli/knob_new_hyd config_vknob_e.xml
arm vknob_f ../cli/knob_old_hyd config_vknob_f.xml
arm vknob_g ../cli/knob_new_hyd config_vknob_g.xml $HK
arm vknob_h ../cli/knob_old_hyd config_vknob_h.xml $HK
wait_arms
for t in a b c d e f g h; do [ -f output_vknob_$t/RUN_CONFIG.txt ] && sed -i 's#output_vknob_[a-h]/#OUT/#' output_vknob_$t/RUN_CONFIG.txt; done
cmp_dirs vknob_a vknob_b "A  atm new clean == old clean"
cmp_dirs vknob_c vknob_d "B  atm new+working branch == old+working branch"
cmp_dirs vknob_e vknob_f "C  hyd new clean == old clean"
cmp_dirs vknob_g vknob_h "D  hyd new+5 knobs == old+5 knobs"
want_differ vknob_c vknob_b "E  CONTROL atm new+working branch vs old clean -- MUST differ"
want_differ vknob_g vknob_f "F  CONTROL hyd new+5 knobs vs old clean -- MUST differ"
echo "--- new banners"; grep -a '\[RUN CONFIG\]' vknob_c.log | cut -c1-250; grep -a '\[RUN CONFIG\]' vknob_g.log | cut -c1-250
touch VKNOB_DONE
