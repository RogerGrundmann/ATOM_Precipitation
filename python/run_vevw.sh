#!/bin/bash
# 2026-10-06: ATM_EVAP_WIND / ATM_EVAP_GUST (defaults 0 / 0 = shipped) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/emlo_atm (7bef4a2 source, -O0; no source change since), new = cli/evwo_atm (HEAD + working tree, -O0). WB = python/working_branch.env
# (wb49 stack; ATM_WATER_CLOSURE=1, so arm B runs the flux branch the knobs act on).
#   A  new clean == old clean      B  new WB == old WB
#   C  CONTROL new WB + ATM_EVAP_WIND=1 != new WB      D  CONTROL new WB + ATM_EVAP_GUST=3 != new WB
# RESULT (10:34): A and B identical 13 of 14 (RUN_CONFIG.txt, banner tokens only); controls C and D differ 13 of 14 -> PASS.
set -u; cd "$(dirname "$0")"; rm -f VEVW_DONE
../python/build_o0.sh WORKTREE evwo || { touch VEVW_DONE; exit 1; }
. ./verify_lib.sh
for t in a c b f d e; do mkdir output_vevw_$t || { touch VEVW_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vevw_$t/#" config_vconv_new.xml > config_vevw_$t.xml; done
arm vevw_a ../cli/evwo_atm config_vevw_a.xml
arm vevw_c ../cli/emlo_atm config_vevw_c.xml
. ./working_branch.env
arm vevw_b ../cli/evwo_atm config_vevw_b.xml
arm vevw_f ../cli/emlo_atm config_vevw_f.xml
arm vevw_d ../cli/evwo_atm config_vevw_d.xml ATM_EVAP_WIND=1
arm vevw_e ../cli/evwo_atm config_vevw_e.xml ATM_EVAP_GUST=3
wait_arms
for t in a c b f d e; do sed -i 's#output_vevw_[a-g]/#OUT/#' output_vevw_$t/RUN_CONFIG.txt; done
cmp_dirs vevw_a vevw_c "A  new clean == old clean"
cmp_dirs vevw_b vevw_f "B  new WB == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  EVAP_(WIND|GUST)=[^ ]*//g" output_vevw_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vevw_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vevw_d vevw_b "C  CONTROL + EVAP_WIND=1 vs off -- MUST differ"
want_differ vevw_e vevw_b "D  CONTROL + EVAP_GUST=3 vs off -- MUST differ"
touch VEVW_DONE
