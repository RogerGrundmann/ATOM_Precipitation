#!/bin/bash
# 2026-10-05: ATM_RH_OCEAN_SST_COLD (default -99 = no cut-off) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/h43o_atm (HEAD 2f4e1bb, -O0), new = cli/ssco_atm (HEAD + working tree, -O0). WB = python/working_branch.env (the wb34 stack).
#   A  new clean == old clean
#   B  new WB + RH_OCEAN_SST=0.007 == old WB + RH_OCEAN_SST=0.007      (the rewritten reduction, cut-off unset)
#   C  CONTROL new WB + RH_OCEAN_SST=0.007 + _COLD=25.5 != new WB + RH_OCEAN_SST=0.007   (the cut-off fires)
# RESULT (12:05): A and B identical 13 of 14 (RUN_CONFIG.txt, banner token only); control C differs 13 of 14 -> PASS.
set -u; cd "$(dirname "$0")"; rm -f VSSC_DONE
. ./verify_lib.sh
for t in a c b f d; do mkdir output_vssc_$t || { touch VSSC_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vssc_$t/#" config_vconv_new.xml > config_vssc_$t.xml; done
arm vssc_a ../cli/ssco_atm config_vssc_a.xml
arm vssc_c ../cli/h43o_atm config_vssc_c.xml
. ./working_branch.env
arm vssc_b ../cli/ssco_atm config_vssc_b.xml ATM_RH_OCEAN_SST=0.007
arm vssc_f ../cli/h43o_atm config_vssc_f.xml ATM_RH_OCEAN_SST=0.007
arm vssc_d ../cli/ssco_atm config_vssc_d.xml ATM_RH_OCEAN_SST=0.007 ATM_RH_OCEAN_SST_COLD=25.5
wait_arms
for t in a c b f d; do sed -i 's#output_vssc_[a-g]/#OUT/#' output_vssc_$t/RUN_CONFIG.txt; done
cmp_dirs vssc_a vssc_c "A  new clean == old clean"
cmp_dirs vssc_b vssc_f "B  new WB RH_OCEAN_SST=0.007 == old WB RH_OCEAN_SST=0.007"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  RH_OCEAN_SST_COLD=[^ ]*//g" output_vssc_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vssc_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vssc_d vssc_b "C  CONTROL + RH_OCEAN_SST_COLD=25.5 vs no cut-off -- MUST differ"
touch VSSC_DONE
