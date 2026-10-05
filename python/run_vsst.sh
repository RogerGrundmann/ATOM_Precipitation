#!/bin/bash
# 2026-10-05: ATM_RH_OCEAN_SST (default 0) + _REF (29.5) + _MAX (0.03) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/h41o_atm (HEAD 972bbb1, -O0), new = cli/ssto_atm (HEAD + working tree, -O0). WB = python/working_branch.env (the wb34 stack).
#   A  new clean == old clean
#   B  new WB + RH_OCEAN_SST=0 == old WB
#   C  CONTROL new WB + RH_OCEAN_SST=0.007 != old WB
#   D  CONTROL new WB + RH_OCEAN_SST=0.007 + _MAX=0.01 != new WB + RH_OCEAN_SST=0.007   (the cap fires)
# RESULT (11:34): A and B identical 13 of 14 (RUN_CONFIG.txt, banner tokens only); controls C and D differ 13 of 14 -> PASS.
set -u; cd "$(dirname "$0")"; rm -f VSST_DONE
. ./verify_lib.sh
for t in a c b f d e; do mkdir output_vsst_$t || { touch VSST_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vsst_$t/#" config_vconv_new.xml > config_vsst_$t.xml; done
arm vsst_a ../cli/ssto_atm config_vsst_a.xml
arm vsst_c ../cli/h41o_atm config_vsst_c.xml
. ./working_branch.env
arm vsst_b ../cli/ssto_atm config_vsst_b.xml ATM_RH_OCEAN_SST=0
arm vsst_f ../cli/h41o_atm config_vsst_f.xml
arm vsst_d ../cli/ssto_atm config_vsst_d.xml ATM_RH_OCEAN_SST=0.007
arm vsst_e ../cli/ssto_atm config_vsst_e.xml ATM_RH_OCEAN_SST=0.007 ATM_RH_OCEAN_SST_MAX=0.01
wait_arms
for t in a c b f d e; do sed -i 's#output_vsst_[a-g]/#OUT/#' output_vsst_$t/RUN_CONFIG.txt; done
cmp_dirs vsst_a vsst_c "A  new clean == old clean"
cmp_dirs vsst_b vsst_f "B  new WB RH_OCEAN_SST=0 == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  RH_OCEAN_SST(_REF|_MAX)?=[^ ]*//g" output_vsst_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vsst_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner tokens" || echo "$1: RUN_CONFIG differs BEYOND the banner tokens"; done
want_differ vsst_d vsst_f "C  CONTROL new WB RH_OCEAN_SST=0.007 vs old WB -- MUST differ"
want_differ vsst_e vsst_d "D  CONTROL + RH_OCEAN_SST_MAX=0.01 vs the default 0.03 -- MUST differ"
touch VSST_DONE
