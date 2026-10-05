#!/bin/bash
# 2026-10-05: ATM_HCRIT_SFC (default 0) + ATM_HCRIT_SFC_LAT (default 0) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/h35o_atm (HEAD 093afd8, -O0), new = cli/hcso_atm (HEAD + working tree, -O0). WB = python/working_branch.env (the wb34 stack).
#   A  new clean == old clean
#   B  new WB + HCRIT_SFC=0 == old WB
#   C  CONTROL new WB + HCRIT_SFC=1 != old WB
#   D  CONTROL new WB + HCRIT_SFC=1 + HCRIT_SFC_LAT=15 != new WB + HCRIT_SFC=1   (the latitude limit fires)
set -u; cd "$(dirname "$0")"; rm -f VHCS_DONE
. ./verify_lib.sh
for t in a c b f d e; do mkdir output_vhcs_$t || { touch VHCS_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vhcs_$t/#" config_vconv_new.xml > config_vhcs_$t.xml; done
arm vhcs_a ../cli/hcso_atm config_vhcs_a.xml
arm vhcs_c ../cli/h35o_atm config_vhcs_c.xml
. ./working_branch.env
arm vhcs_b ../cli/hcso_atm config_vhcs_b.xml ATM_HCRIT_SFC=0
arm vhcs_f ../cli/h35o_atm config_vhcs_f.xml
arm vhcs_d ../cli/hcso_atm config_vhcs_d.xml ATM_HCRIT_SFC=1
arm vhcs_e ../cli/hcso_atm config_vhcs_e.xml ATM_HCRIT_SFC=1 ATM_HCRIT_SFC_LAT=15
wait_arms
for t in a c b f d e; do sed -i 's#output_vhcs_[a-g]/#OUT/#' output_vhcs_$t/RUN_CONFIG.txt; done
cmp_dirs vhcs_a vhcs_c "A  new clean == old clean"
cmp_dirs vhcs_b vhcs_f "B  new WB HCRIT_SFC=0 == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  HCRIT_SFC(_LAT)?=[^ ]*//g" output_vhcs_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vhcs_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner tokens" || echo "$1: RUN_CONFIG differs BEYOND the banner tokens"; done
want_differ vhcs_d vhcs_f "C  CONTROL new WB HCRIT_SFC=1 vs old WB -- MUST differ"
want_differ vhcs_e vhcs_d "D  CONTROL + HCRIT_SFC_LAT=15 vs every latitude -- MUST differ"
touch VHCS_DONE
