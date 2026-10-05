#!/bin/bash
# 2026-10-05: ATM_HCRIT_SFC + ATM_HCRIT_SFC_LAT, SECOND form (ground >= 1000 hPa left alone) -- -O0 byte check, 1 thread, nm 20 from scratch. run_vhcs.sh checked the first form (passed 13 of 14, banner tokens only).
# RESULT (11:10): A and B identical 13 of 14 (RUN_CONFIG.txt, banner tokens only); controls C and D differ 13 of 14 -> PASS.
# old = cli/h35o_atm (HEAD 093afd8, -O0), new = cli/hc2o_atm (HEAD + working tree, -O0). WB = python/working_branch.env (the wb34 stack).
#   A  new clean == old clean
#   B  new WB + HCRIT_SFC=0 == old WB
#   C  CONTROL new WB + HCRIT_SFC=1 != old WB
#   D  CONTROL new WB + HCRIT_SFC=1 + HCRIT_SFC_LAT=15 != new WB + HCRIT_SFC=1   (the latitude limit fires)
set -u; cd "$(dirname "$0")"; rm -f VHC2_DONE
. ./verify_lib.sh
for t in a c b f d e; do mkdir output_vhc2_$t || { touch VHC2_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vhc2_$t/#" config_vconv_new.xml > config_vhc2_$t.xml; done
arm vhc2_a ../cli/hc2o_atm config_vhc2_a.xml
arm vhc2_c ../cli/h35o_atm config_vhc2_c.xml
. ./working_branch.env
arm vhc2_b ../cli/hc2o_atm config_vhc2_b.xml ATM_HCRIT_SFC=0
arm vhc2_f ../cli/h35o_atm config_vhc2_f.xml
arm vhc2_d ../cli/hc2o_atm config_vhc2_d.xml ATM_HCRIT_SFC=1
arm vhc2_e ../cli/hc2o_atm config_vhc2_e.xml ATM_HCRIT_SFC=1 ATM_HCRIT_SFC_LAT=15
wait_arms
for t in a c b f d e; do sed -i 's#output_vhc2_[a-g]/#OUT/#' output_vhc2_$t/RUN_CONFIG.txt; done
cmp_dirs vhc2_a vhc2_c "A  new clean == old clean"
cmp_dirs vhc2_b vhc2_f "B  new WB HCRIT_SFC=0 == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  HCRIT_SFC(_LAT)?=[^ ]*//g" output_vhc2_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vhc2_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner tokens" || echo "$1: RUN_CONFIG differs BEYOND the banner tokens"; done
want_differ vhc2_d vhc2_f "C  CONTROL new WB HCRIT_SFC=1 vs old WB -- MUST differ"
want_differ vhc2_e vhc2_d "D  CONTROL + HCRIT_SFC_LAT=15 vs every latitude -- MUST differ"
touch VHC2_DONE
