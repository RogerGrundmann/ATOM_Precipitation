#!/bin/bash
# 2026-10-02: ATM_RH_OCEAN (default 0.75 = shipped) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/blqo_atm (364af9d), new = cli/rhoo_atm (+ the knob).
#   A  new clean == old clean;  B  new MC_knobs=0 == old clean;  C  CONTROL new MC_knobs=1 != old clean
set -u; cd "$(dirname "$0")"; rm -f VRHO_DONE
. ./verify_lib.sh
for t in a b c d; do mkdir output_vrho_$t || { touch VRHO_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vrho_$t/#" config_vconv_new.xml > config_vrho_$t.xml; done
arm vrho_a ../cli/rhoo_atm config_vrho_a.xml
arm vrho_b ../cli/rhoo_atm config_vrho_b.xml ATM_RH_OCEAN=0.75
arm vrho_c ../cli/blqo_atm  config_vrho_c.xml
arm vrho_d ../cli/rhoo_atm config_vrho_d.xml ATM_RH_OCEAN=0.82
wait_arms
for t in a b c d; do sed -i 's#output_vrho_[a-d]/#OUT/#' output_vrho_$t/RUN_CONFIG.txt; done
cmp_dirs vrho_a vrho_c "A  new clean == old clean"
cmp_dirs vrho_b vrho_c "B  new MC_knobs=0 == old clean"
for t in a b; do diff <(sed -E "s/  RH_OCEAN=[^ ]*//g" output_vrho_$t/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vrho_c/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$t: RUN_CONFIG differs only by the banner token" || echo "$t: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vrho_d vrho_c "C  CONTROL new MC_knobs=1 vs old clean -- MUST differ"
touch VRHO_DONE
