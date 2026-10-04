#!/bin/bash
# 2026-10-04: ATM_RH_OCEAN_ML (default 0) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/h28o_atm (HEAD 211fd08, -O0), new = cli/omlo_atm (HEAD + working tree, -O0). WB = python/working_branch.env.
#   A  new clean == old clean
#   B  new WB + RH_OCEAN_ML=0 == old WB
#   C  CONTROL new WB + RH_OCEAN_ML=1200 != old WB
set -u; cd "$(dirname "$0")"; rm -f VMOML_DONE
. ./verify_lib.sh
for t in a c b f d; do mkdir output_vmoml_$t || { touch VMOML_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vmoml_$t/#" config_vconv_new.xml > config_vmoml_$t.xml; done
arm vmoml_a ../cli/omlo_atm  config_vmoml_a.xml
arm vmoml_c ../cli/h28o_atm config_vmoml_c.xml
. ./working_branch.env
arm vmoml_b ../cli/omlo_atm  config_vmoml_b.xml ATM_RH_OCEAN_ML=0
arm vmoml_f ../cli/h28o_atm config_vmoml_f.xml
arm vmoml_d ../cli/omlo_atm  config_vmoml_d.xml ATM_RH_OCEAN_ML=1200
wait_arms
for t in a c b f d; do sed -i 's#output_vmoml_[a-g]/#OUT/#' output_vmoml_$t/RUN_CONFIG.txt; done
grep -a "nm\b\|nm=" output_vmoml_a/RUN_CONFIG.txt | head -2
cmp_dirs vmoml_a vmoml_c "A  new clean == old clean"
cmp_dirs vmoml_b vmoml_f "B  new WB RH_OCEAN_ML=0 == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  RH_OCEAN_ML=[^ ]*//g" output_vmoml_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vmoml_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vmoml_d vmoml_f "C  CONTROL new WB RH_OCEAN_ML=1200 vs old WB -- MUST differ"
touch VMOML_DONE

