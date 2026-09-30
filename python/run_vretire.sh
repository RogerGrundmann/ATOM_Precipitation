#!/bin/bash
# KNOB-INV plan C: retiring decided knobs must change NOTHING on the default branch.
# usage: run_vretire.sh <tag> <old binary prefix, e.g. knob_new>   (new = HEAD + working tree, -O0)
# -O0 both sides, 1 thread.  atm nm 20 from scratch: A clean, B working_branch.env;  hyd 1000 -> 1020 from oc_ctl: C clean.
# Controls: E atm new+working branch vs old clean MUST differ.  RUN_CONFIG.txt may differ (banner shrinks).
set -u; cd "$(dirname "$0")"; T=$1; OLD=$2; rm -f VRET_${T}_DONE
. ./verify_lib.sh
./build_o0.sh WORKTREE ret_$T || { touch VRET_${T}_DONE; exit 1; }
WB=$(grep '^export ' working_branch.env | sed 's/^export //; s/ *#.*$//' | tr '\n' ' ')
for t in a b c d; do mkdir output_vret${T}_$t || { touch VRET_${T}_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vret${T}_$t/#" config_vconv_new.xml > config_vret${T}_$t.xml; done
for t in e f; do mkdir output_vret${T}_$t || { touch VRET_${T}_DONE; exit 1; }
  sed -e "s#output_sm_bc1/#output_vret${T}_$t/#" -e "s#<nm>200</nm>#<nm>20</nm>#" \
      -e "s#<checkpoint_save_iter>1200</checkpoint_save_iter>#<checkpoint_save_iter>1020</checkpoint_save_iter>#" config_sm_bc1.xml > config_vret${T}_$t.xml
  cp output_oc_ctl/hyd_restart_0Ma_1000.bin output_vret${T}_$t/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_vret${T}_$t/; done
arm vret${T}_a ../cli/ret_${T}_atm config_vret${T}_a.xml
arm vret${T}_b ../cli/${OLD}_atm   config_vret${T}_b.xml
arm vret${T}_c ../cli/ret_${T}_atm config_vret${T}_c.xml $WB
arm vret${T}_d ../cli/${OLD}_atm   config_vret${T}_d.xml $WB
arm vret${T}_e ../cli/ret_${T}_hyd config_vret${T}_e.xml
arm vret${T}_f ../cli/${OLD}_hyd   config_vret${T}_f.xml
wait_arms
for t in a b c d e f; do [ -f output_vret${T}_$t/RUN_CONFIG.txt ] && sed -i "s#output_vret${T}_[a-f]/#OUT/#" output_vret${T}_$t/RUN_CONFIG.txt; done
cmp_dirs vret${T}_a vret${T}_b "A  atm new clean == old clean"
cmp_dirs vret${T}_c vret${T}_d "B  atm new+working branch == old+working branch"
cmp_dirs vret${T}_e vret${T}_f "C  hyd new clean == old clean"
want_differ vret${T}_c vret${T}_b "E  CONTROL atm new+working branch vs old clean -- MUST differ"
BAD=$(grep 'DIFFERS:' run_vretire_$T.out | grep -vc 'RUN_CONFIG.txt')
echo "VERDICT: $([ "$BAD" = 0 ] && grep -q 'E  CONTROL.*PASS' run_vretire_$T.out && echo PASS || echo FAIL)"
touch VRET_${T}_DONE
