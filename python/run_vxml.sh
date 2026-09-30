#!/bin/bash
# 2026-09-30: KNOB-INV plan D -- knobs from a <knobs> section of the XML config (env > XML > default).
# -O0 both sides, 1 thread. old = cli/ret_b5_* (bc4e63a), new = HEAD + working tree.
#   A atm clean == old clean;  B atm +working_branch env == old +env;  X atm with the SAME values in <knobs> == old +env
#   C hyd clean == old clean;  Y hyd with 5 knobs in <knobs> == old with them in the env
#   E CONTROL atm new+env vs old clean MUST differ;  T a <knobs> typo must make the run FAIL (non-zero exit)
set -u; cd "$(dirname "$0")"; rm -f VXML_DONE
. ./verify_lib.sh
./build_o0.sh WORKTREE xml || { touch VXML_DONE; exit 1; }
WB=$(grep '^export ' working_branch.env | sed 's/^export //; s/ *#.*$//' | tr '\n' ' ')
HK="HYD_METRIC_RADIUS=0 HYD_RUN_NEUMANN=0 HYD_A_H_BIHARM=0 HYD_A_H_BIHARM_SCALED=0 HYD_T_FREEZE_SFC=1"
knobs_xml(){ echo "  <knobs>"; for kv in $1; do echo "    <${kv%%=*}>${kv#*=}</${kv%%=*}>"; done; echo "  </knobs>"; }
for t in a b c d x t; do mkdir output_vxml_$t || { touch VXML_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vxml_$t/#" config_vconv_new.xml > config_vxml_$t.xml; done
python3 - "$(knobs_xml "$WB")" config_vxml_x.xml <<'PY'
import sys; blk, p = sys.argv[1], sys.argv[2]; s = open(p).read(); assert s.count('</atom>') == 1
open(p, 'w').write(s.replace('</atom>', blk + '\n</atom>'))
PY
python3 - "$(knobs_xml "ATM_RH_STORM=1.15 ATM_NOT_A_KNOB=1")" config_vxml_t.xml <<'PY'
import sys; blk, p = sys.argv[1], sys.argv[2]; s = open(p).read(); open(p, 'w').write(s.replace('</atom>', blk + '\n</atom>'))
PY
for t in e f y z; do mkdir output_vxml_$t || { touch VXML_DONE; exit 1; }
  sed -e "s#output_sm_bc1/#output_vxml_$t/#" -e "s#<nm>200</nm>#<nm>20</nm>#" \
      -e "s#<checkpoint_save_iter>1200</checkpoint_save_iter>#<checkpoint_save_iter>1020</checkpoint_save_iter>#" config_sm_bc1.xml > config_vxml_$t.xml
  cp output_oc_ctl/hyd_restart_0Ma_1000.bin output_vxml_$t/; cp output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp output_vxml_$t/; done
python3 - "$(knobs_xml "$HK")" config_vxml_y.xml <<'PY'
import sys; blk, p = sys.argv[1], sys.argv[2]; s = open(p).read(); assert s.count('</atom>') == 1
open(p, 'w').write(s.replace('</atom>', blk + '\n</atom>'))
PY
arm vxml_a ../cli/xml_atm     config_vxml_a.xml
arm vxml_b ../cli/ret_b5_atm  config_vxml_b.xml
arm vxml_c ../cli/xml_atm     config_vxml_c.xml $WB
arm vxml_d ../cli/ret_b5_atm  config_vxml_d.xml $WB
arm vxml_x ../cli/xml_atm     config_vxml_x.xml
arm vxml_e ../cli/xml_hyd     config_vxml_e.xml
arm vxml_f ../cli/ret_b5_hyd  config_vxml_f.xml
arm vxml_y ../cli/xml_hyd     config_vxml_y.xml
arm vxml_z ../cli/ret_b5_hyd  config_vxml_z.xml $HK
wait_arms
env OMP_NUM_THREADS=1 ../cli/xml_atm config_vxml_t.xml > vxml_t.log 2>&1; te=$?
for t in a b c d x e f y z; do [ -f output_vxml_$t/RUN_CONFIG.txt ] && sed -i "s#output_vxml_[a-z]/#OUT/#" output_vxml_$t/RUN_CONFIG.txt; done
cmp_dirs vxml_a vxml_b "A  atm new clean == old clean"
cmp_dirs vxml_c vxml_d "B  atm new+env == old+env"
cmp_dirs vxml_x vxml_d "X  atm new with <knobs> == old+env"
cmp_dirs vxml_e vxml_f "C  hyd new clean == old clean"
cmp_dirs vxml_y vxml_z "D  hyd new with <knobs> == old+env"
want_differ vxml_c vxml_b "E  CONTROL atm new+env vs old clean -- MUST differ"
echo "T  typo in <knobs>: exit $te -> $([ $te != 0 ] && echo PASS || echo FAIL); $(grep -a -m1 'unknown knob' vxml_t.log | cut -c1-140)"
echo "--- banner of the XML arm:"; grep -a '\[RUN CONFIG\]' vxml_x.log | grep -a 'RH_STORM\|+ = from' | cut -c1-200
touch VXML_DONE
