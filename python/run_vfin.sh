#!/bin/bash
# 2026-10-03: FINAL-TREE re-check after run_vct.sh failed B at round-off level (v-budget terms ~1e-19) with the distance Dijkstra MOVED out of
# initWaterWapour. The block is back in place (LandDistance.h is a copy for MoistConvection). old = cli/tlo_atm (88ad872 + ATM_MC_T_ADD_LAND),
# new = cli/fin_atm (+ ATM_MC_T_ADD_LAND_L + ATM_MC_MB_SAT_LAND). K = ATM_MC_ML_PARCEL=1 ATM_RH_LAND=2 ATM_MC_T_ADD_LAND=2.
#   A  new clean == old clean;  B  new K + L=0 + MB_SAT_LAND=0 == old K;  C  CONTROL new K + L=300 != old K
set -u; cd "$(dirname "$0")"; rm -f VFIN_DONE
. ./verify_lib.sh
for t in a b c d f; do mkdir output_vfin_$t || { touch VFIN_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vfin_$t/#" config_vconv_new.xml > config_vfin_$t.xml; done
K="ATM_MC_ML_PARCEL=1 ATM_RH_LAND=2 ATM_MC_T_ADD_LAND=2"
arm vfin_a ../cli/fin_atm config_vfin_a.xml
arm vfin_c ../cli/tlo_atm config_vfin_c.xml
arm vfin_b ../cli/fin_atm config_vfin_b.xml $K ATM_MC_T_ADD_LAND_L=0 ATM_MC_MB_SAT_LAND=0
arm vfin_f ../cli/tlo_atm config_vfin_f.xml $K
arm vfin_d ../cli/fin_atm config_vfin_d.xml $K ATM_MC_T_ADD_LAND_L=300
wait_arms
for t in a b c d f; do sed -i 's#output_vfin_[a-f]/#OUT/#' output_vfin_$t/RUN_CONFIG.txt; done
cmp_dirs vfin_a vfin_c "A  new clean == old clean"
cmp_dirs vfin_b vfin_f "B  new K L=0 == old K"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  MC_(T_ADD_LAND_L|MB_SAT_LAND)=[^ ]*//g" output_vfin_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vfin_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vfin_d vfin_f "C  CONTROL new K L=300 vs old K -- MUST differ"
touch VFIN_DONE
