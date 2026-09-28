#!/bin/bash
# 2026-09-28: both-directions -O0 byte check for the ATM_TURB_SIN_FLOOR flip (branch tsf-flip, 1190b14).
# old = cli/conv_new_atm (c1aae23 -O0), new = cli/tsf_new_atm (tsf-flip -O0). 1 thread, nm 20 from scratch.
#   A  new clean            == old ATM_TURB_SIN_FLOOR=1   (the flip reproduces the knob)
#   B  new ATM_TURB_SIN_FLOOR=0 == old clean               (the old branch is restorable)
#   C  new clean            != old clean                   (control: the flip fires at 20 iterations)
# RUN_CONFIG.txt may differ only by the banner token (1* vs 1 / 0 vs 0*); printed. Starts after run_rp0928.sh.
set -u; cd "$(dirname "$0")"; rm -f VTSF0928_DONE
until [ -e RP0928_DONE ]; do sleep 60; done; sleep 240   # stagger: run_sqd0928.sh builds its -O0 worktree at the same moment
. ./verify_lib.sh
[ -e ../cli/tsf_new_atm ] || ./build_o0.sh tsf-flip tsf_new || { touch VTSF0928_DONE; exit 1; }
for t in a b c d; do mkdir output_vtsf_$t || { touch VTSF0928_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vtsf_$t/#" config_vconv_new.xml > config_vtsf_$t.xml; done
arm vtsf_a ../cli/tsf_new_atm  config_vtsf_a.xml
arm vtsf_b ../cli/conv_new_atm config_vtsf_b.xml ATM_TURB_SIN_FLOOR=1
arm vtsf_c ../cli/tsf_new_atm  config_vtsf_c.xml ATM_TURB_SIN_FLOOR=0
arm vtsf_d ../cli/conv_new_atm config_vtsf_d.xml
wait_arms
for t in a b c d; do sed -i 's#output_vtsf_[a-d]/#OUT/#' output_vtsf_$t/RUN_CONFIG.txt; done
cmp_dirs vtsf_a vtsf_b "A  FLIP == old+knob"
diff output_vtsf_a/RUN_CONFIG.txt output_vtsf_b/RUN_CONFIG.txt | grep -o 'TURB_SIN_FLOOR=[^ ]*'
cmp_dirs vtsf_c vtsf_d "B  new+knob=0 == old clean"
diff output_vtsf_c/RUN_CONFIG.txt output_vtsf_d/RUN_CONFIG.txt | grep -o 'TURB_SIN_FLOOR=[^ ]*'
want_differ vtsf_a vtsf_d "C  CONTROL new clean vs old clean -- MUST differ"
touch VTSF0928_DONE
