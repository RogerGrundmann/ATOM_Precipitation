#!/bin/bash
# 2026-10-02: ATM_RH_LAND (default 0) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/lclc_atm (20319f8), new = cli/rhlo_atm (+ the knob).
#   A  new clean == old clean;  B  new MC_knobs=0 == old clean;  C  CONTROL new MC_knobs=1 != old clean
set -u; cd "$(dirname "$0")"; rm -f VRHL_DONE
. ./verify_lib.sh
for t in a b c d; do mkdir output_vrhl_$t || { touch VRHL_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vrhl_$t/#" config_vconv_new.xml > config_vrhl_$t.xml; done
arm vrhl_a ../cli/rhlo_atm config_vrhl_a.xml
arm vrhl_b ../cli/rhlo_atm config_vrhl_b.xml ATM_RH_LAND=0
arm vrhl_c ../cli/lclc_atm  config_vrhl_c.xml
arm vrhl_d ../cli/rhlo_atm config_vrhl_d.xml ATM_RH_LAND=1
wait_arms
for t in a b c d; do sed -i 's#output_vrhl_[a-d]/#OUT/#' output_vrhl_$t/RUN_CONFIG.txt; done
cmp_dirs vrhl_a vrhl_c "A  new clean == old clean"
cmp_dirs vrhl_b vrhl_c "B  new MC_knobs=0 == old clean"
for t in a b; do diff <(sed -E "s/  RH_LAND=[^ ]*//g" output_vrhl_$t/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vrhl_c/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$t: RUN_CONFIG differs only by the banner token" || echo "$t: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vrhl_d vrhl_c "C  CONTROL new MC_knobs=1 vs old clean -- MUST differ"
touch VRHL_DONE
