#!/bin/bash
# 2026-10-04: ATM_RH_OCEAN_ML_STRENGTH (default 1) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/omlo_atm (= d8193cd source, -O0), new = cli/omso_atm (d8193cd + working tree, -O0). WB = python/working_branch.env.
#   A  new clean == old clean
#   B  new WB + RH_OCEAN_ML=1200 (strength at its default 1) == old WB + RH_OCEAN_ML=1200      (the ON branch is unchanged)
#   C  CONTROL new WB + RH_OCEAN_ML=1200 + STRENGTH=0.4 != old WB + RH_OCEAN_ML=1200
set -u; cd "$(dirname "$0")"; rm -f VMOMS_DONE
. ./verify_lib.sh
for t in a c b f d; do mkdir output_vmoms_$t || { touch VMOMS_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vmoms_$t/#" config_vconv_new.xml > config_vmoms_$t.xml; done
arm vmoms_a ../cli/omso_atm config_vmoms_a.xml
arm vmoms_c ../cli/omlo_atm config_vmoms_c.xml
. ./working_branch.env
arm vmoms_b ../cli/omso_atm config_vmoms_b.xml ATM_RH_OCEAN_ML=1200
arm vmoms_f ../cli/omlo_atm config_vmoms_f.xml ATM_RH_OCEAN_ML=1200
arm vmoms_d ../cli/omso_atm config_vmoms_d.xml ATM_RH_OCEAN_ML=1200 ATM_RH_OCEAN_ML_STRENGTH=0.4
wait_arms
for t in a c b f d; do sed -i 's#output_vmoms_[a-g]/#OUT/#' output_vmoms_$t/RUN_CONFIG.txt; done
cmp_dirs vmoms_a vmoms_c "A  new clean == old clean"
cmp_dirs vmoms_b vmoms_f "B  new WB ML=1200 == old WB ML=1200"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  RH_OCEAN_ML_STRENGTH=[^ ]*//g" output_vmoms_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vmoms_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vmoms_d vmoms_f "C  CONTROL new WB ML=1200 STRENGTH=0.4 vs old WB ML=1200 -- MUST differ"
touch VMOMS_DONE
