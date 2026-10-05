#!/bin/bash
# 2026-10-05: ATM_RH_LAND_EAST_MAX (default -1 = no clamp) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/h47o_atm (HEAD 12ed221, -O0), new = cli/emxo_atm (HEAD + working tree, -O0). WB = python/working_branch.env (the wb45 stack,
# which has ATM_RH_LAND_EAST=1.35 -- so arm B exercises the rewritten return of rh_land_of with the clamp unset).
#   A  new clean == old clean
#   B  new WB == old WB
#   C  CONTROL new WB + RH_LAND_EAST_MAX=0.015 != new WB   (the clamp fires)
# RESULT (14:38): A and B identical 13 of 14 (RUN_CONFIG.txt, banner token only); control C differs 12 of 14 -> PASS.
set -u; cd "$(dirname "$0")"; rm -f VEMX_DONE
. ./verify_lib.sh
for t in a c b f d; do mkdir output_vemx_$t || { touch VEMX_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vemx_$t/#" config_vconv_new.xml > config_vemx_$t.xml; done
arm vemx_a ../cli/emxo_atm config_vemx_a.xml
arm vemx_c ../cli/h47o_atm config_vemx_c.xml
. ./working_branch.env
arm vemx_b ../cli/emxo_atm config_vemx_b.xml
arm vemx_f ../cli/h47o_atm config_vemx_f.xml
arm vemx_d ../cli/emxo_atm config_vemx_d.xml ATM_RH_LAND_EAST_MAX=0.015
wait_arms
for t in a c b f d; do sed -i 's#output_vemx_[a-g]/#OUT/#' output_vemx_$t/RUN_CONFIG.txt; done
cmp_dirs vemx_a vemx_c "A  new clean == old clean"
cmp_dirs vemx_b vemx_f "B  new WB == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  RH_LAND_EAST_MAX=[^ ]*//g" output_vemx_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vemx_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vemx_d vemx_b "C  CONTROL + RH_LAND_EAST_MAX=0.015 vs no clamp -- MUST differ"
touch VEMX_DONE
