#!/bin/bash
# 2026-10-03: ATM_MC_MB_SAT_LAND (default 0) -- -O0 byte check, 1 thread, nm 20 from scratch. Waits for the vct check (its new side is this old side).
# old = cli/cto_atm (88ad872 + ATM_MC_T_ADD_LAND + _L), new = cli/mbo_atm (+ the knob). K = ATM_MC_ML_PARCEL=1 (land convection active).
#   A  new clean == old clean;  B  new K + MB_SAT_LAND=0 == old K;  C  CONTROL new K + MB_SAT_LAND=0.04 != old K
set -u; cd "$(dirname "$0")"; rm -f VMB_DONE
until [ -e VCT_DONE ]; do sleep 20; done
. ./verify_lib.sh
for t in a b c d f; do mkdir output_vmb_$t || { touch VMB_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vmb_$t/#" config_vconv_new.xml > config_vmb_$t.xml; done
K="ATM_MC_ML_PARCEL=1"
arm vmb_a ../cli/mbo_atm config_vmb_a.xml
arm vmb_c ../cli/cto_atm config_vmb_c.xml
arm vmb_b ../cli/mbo_atm config_vmb_b.xml $K ATM_MC_MB_SAT_LAND=0
arm vmb_f ../cli/cto_atm config_vmb_f.xml $K
arm vmb_d ../cli/mbo_atm config_vmb_d.xml $K ATM_MC_MB_SAT_LAND=0.04
wait_arms
for t in a b c d f; do sed -i 's#output_vmb_[a-f]/#OUT/#' output_vmb_$t/RUN_CONFIG.txt; done
cmp_dirs vmb_a vmb_c "A  new clean == old clean"
cmp_dirs vmb_b vmb_f "B  new K MB_SAT_LAND=0 == old K"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  MC_MB_SAT_LAND=[^ ]*//g" output_vmb_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vmb_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vmb_d vmb_f "C  CONTROL new K MB_SAT_LAND=0.04 vs old K -- MUST differ"
touch VMB_DONE
