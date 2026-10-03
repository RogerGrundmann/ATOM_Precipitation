#!/bin/bash
# 2026-10-03: ATM_MC_T_ADD_LAND_L (default 0) + the distance-to-ocean Dijkstra moved to LandDistance.h -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/tlo_atm (88ad872 + ATM_MC_T_ADD_LAND), new = cli/cto_atm. K = ATM_MC_ML_PARCEL=1 ATM_RH_LAND=2 ATM_MC_T_ADD_LAND=2 exercises both the
# moved Dijkstra (ATM_RH_LAND path) and the knob's off branch.
#   A  new clean == old clean;  B  new K + L=0 == old K;  C  CONTROL new K + L=300 != old K
set -u; cd "$(dirname "$0")"; rm -f VCT_DONE
. ./verify_lib.sh
for t in a b c d f; do mkdir output_vct_$t || { touch VCT_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vct_$t/#" config_vconv_new.xml > config_vct_$t.xml; done
K="ATM_MC_ML_PARCEL=1 ATM_RH_LAND=2 ATM_MC_T_ADD_LAND=2"
arm vct_a ../cli/cto_atm config_vct_a.xml
arm vct_c ../cli/tlo_atm config_vct_c.xml
arm vct_b ../cli/cto_atm config_vct_b.xml $K ATM_MC_T_ADD_LAND_L=0
arm vct_f ../cli/tlo_atm config_vct_f.xml $K
arm vct_d ../cli/cto_atm config_vct_d.xml $K ATM_MC_T_ADD_LAND_L=300
wait_arms
for t in a b c d f; do sed -i 's#output_vct_[a-f]/#OUT/#' output_vct_$t/RUN_CONFIG.txt; done
cmp_dirs vct_a vct_c "A  new clean == old clean"
cmp_dirs vct_b vct_f "B  new K L=0 == old K"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  MC_T_ADD_LAND_L=[^ ]*//g" output_vct_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vct_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vct_d vct_f "C  CONTROL new K L=300 vs old K -- MUST differ"
touch VCT_DONE
