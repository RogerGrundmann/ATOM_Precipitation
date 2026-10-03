#!/bin/bash
# 2026-10-03: ATM_MC_GATE_BLEND (default 0) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/gbh_atm (4960020), new = cli/gbo_atm (+ the knob).
#   A  new clean == old clean;  B  new knob=0 == old clean;  C  CONTROL new knob=1 != old clean;  D  CONTROL knob=2 != knob=1
set -u; cd "$(dirname "$0")"; rm -f VGB_DONE
. ./verify_lib.sh
for t in a b c d e; do mkdir output_vgb_$t || { touch VGB_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vgb_$t/#" config_vconv_new.xml > config_vgb_$t.xml; done
arm vgb_a ../cli/gbo_atm config_vgb_a.xml
arm vgb_b ../cli/gbo_atm config_vgb_b.xml ATM_MC_GATE_BLEND=0
arm vgb_c ../cli/gbh_atm config_vgb_c.xml
arm vgb_d ../cli/gbo_atm config_vgb_d.xml ATM_MC_GATE_BLEND=1
arm vgb_e ../cli/gbo_atm config_vgb_e.xml ATM_MC_GATE_BLEND=2
wait_arms
for t in a b c d e; do sed -i 's#output_vgb_[a-e]/#OUT/#' output_vgb_$t/RUN_CONFIG.txt; done
cmp_dirs vgb_a vgb_c "A  new clean == old clean"
cmp_dirs vgb_b vgb_c "B  new knob=0 == old clean"
for t in a b; do diff <(sed -E "s/  MC_GATE_BLEND=[^ ]*//g" output_vgb_$t/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vgb_c/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$t: RUN_CONFIG differs only by the banner token" || echo "$t: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vgb_d vgb_c "C  CONTROL new knob=1 vs old clean -- MUST differ"
want_differ vgb_e vgb_d "D  CONTROL knob=2 vs knob=1 -- MUST differ"
touch VGB_DONE
