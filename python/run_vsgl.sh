#!/bin/bash
# 2026-10-04: ATM_RH_SIGMA_LAT (default 0) + ATM_RH_SIGMA_LAT_OCEAN (default 0) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/h34o_atm (HEAD 123e4f0, -O0), new = cli/sglo_atm (HEAD + working tree, -O0). WB = python/working_branch.env (the wb34 stack).
#   A  new clean == old clean
#   B  new WB + RH_SIGMA_LAT=0 == old WB
#   C  CONTROL new WB + RH_SIGMA_LAT=52 != old WB
#   D  CONTROL new WB + RH_SIGMA_LAT=52 + RH_SIGMA_LAT_OCEAN=1 != new WB + RH_SIGMA_LAT=52   (the ocean switch fires)
set -u; cd "$(dirname "$0")"; rm -f VSGL_DONE
. ./verify_lib.sh
for t in a c b f d e; do mkdir output_vsgl_$t || { touch VSGL_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vsgl_$t/#" config_vconv_new.xml > config_vsgl_$t.xml; done
arm vsgl_a ../cli/sglo_atm config_vsgl_a.xml
arm vsgl_c ../cli/h34o_atm config_vsgl_c.xml
. ./working_branch.env
arm vsgl_b ../cli/sglo_atm config_vsgl_b.xml ATM_RH_SIGMA_LAT=0
arm vsgl_f ../cli/h34o_atm config_vsgl_f.xml
arm vsgl_d ../cli/sglo_atm config_vsgl_d.xml ATM_RH_SIGMA_LAT=52
arm vsgl_e ../cli/sglo_atm config_vsgl_e.xml ATM_RH_SIGMA_LAT=52 ATM_RH_SIGMA_LAT_OCEAN=1
wait_arms
for t in a c b f d e; do sed -i 's#output_vsgl_[a-g]/#OUT/#' output_vsgl_$t/RUN_CONFIG.txt; done
cmp_dirs vsgl_a vsgl_c "A  new clean == old clean"
cmp_dirs vsgl_b vsgl_f "B  new WB RH_SIGMA_LAT=0 == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  RH_SIGMA_LAT(_OCEAN)?=[^ ]*//g" output_vsgl_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vsgl_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner tokens" || echo "$1: RUN_CONFIG differs BEYOND the banner tokens"; done
want_differ vsgl_d vsgl_f "C  CONTROL new WB RH_SIGMA_LAT=52 vs old WB -- MUST differ"
want_differ vsgl_e vsgl_d "D  CONTROL + RH_SIGMA_LAT_OCEAN=1 vs land only -- MUST differ"
touch VSGL_DONE
