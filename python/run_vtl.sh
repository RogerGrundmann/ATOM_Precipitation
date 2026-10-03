#!/bin/bash
# 2026-10-03: ATM_MC_T_ADD_LAND (default 0) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/eao_atm (ad25105 source), new = cli/tlo_atm (+ the knob). The knob lives in the ATM_MC_ML_PARCEL path, which the default never
# runs, so the off branch is checked on the default (A) and with ATM_MC_ML_PARCEL=1 (B).
#   A  new clean == old clean;  B  new ML_PARCEL=1 T_ADD_LAND=0 == old ML_PARCEL=1;  C  CONTROL new ML_PARCEL=1 T_ADD_LAND=2 != old ML_PARCEL=1
set -u; cd "$(dirname "$0")"; rm -f VTL_DONE
. ./verify_lib.sh
for t in a b c d f; do mkdir output_vtl_$t || { touch VTL_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vtl_$t/#" config_vconv_new.xml > config_vtl_$t.xml; done
arm vtl_a ../cli/tlo_atm config_vtl_a.xml
arm vtl_c ../cli/eao_atm config_vtl_c.xml
arm vtl_b ../cli/tlo_atm config_vtl_b.xml ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD_LAND=0
arm vtl_f ../cli/eao_atm config_vtl_f.xml ATM_MC_ML_PARCEL=1
arm vtl_d ../cli/tlo_atm config_vtl_d.xml ATM_MC_ML_PARCEL=1 ATM_MC_T_ADD_LAND=2
wait_arms
for t in a b c d f; do sed -i 's#output_vtl_[a-f]/#OUT/#' output_vtl_$t/RUN_CONFIG.txt; done
cmp_dirs vtl_a vtl_c "A  new clean == old clean"
cmp_dirs vtl_b vtl_f "B  new ML_PARCEL=1 T_ADD_LAND=0 == old ML_PARCEL=1"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  MC_T_ADD_LAND=[^ ]*//g" output_vtl_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vtl_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vtl_d vtl_f "C  CONTROL new ML_PARCEL=1 T_ADD_LAND=2 vs old ML_PARCEL=1 -- MUST differ"
touch VTL_DONE
