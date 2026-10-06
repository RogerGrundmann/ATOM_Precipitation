#!/bin/bash
# 2026-10-06: ATM_RH_STORM_LAT / ATM_RH_STORM_WIDTH (defaults 55 / 15 = shipped) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/levo_atm (HEAD + ATM_LAND_EVAP, -O0), new = cli/stoo_atm (the same + the two knobs, -O0). WB = python/working_branch.env (wb57 stack,
# ATM_RH_STORM=1.15 and ATM_RH_STORM_POLAR=1.08, so arm B runs both rewritten expressions).
#   A  new clean == old clean      B  new WB == old WB      C  CONTROL new WB + ATM_RH_STORM_LAT=48 != new WB
# RESULT (12:42): A and B identical 13 of 14 (RUN_CONFIG.txt, banner tokens only); control C differs 13 of 14 -> PASS.
set -u; cd "$(dirname "$0")"; rm -f VSTO_DONE
../python/build_o0.sh WORKTREE stoo || { touch VSTO_DONE; exit 1; }
. ./verify_lib.sh
for t in a c b f d; do mkdir output_vsto_$t || { touch VSTO_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vsto_$t/#" config_vconv_new.xml > config_vsto_$t.xml; done
arm vsto_a ../cli/stoo_atm config_vsto_a.xml
arm vsto_c ../cli/levo_atm config_vsto_c.xml
. ./working_branch.env
arm vsto_b ../cli/stoo_atm config_vsto_b.xml
arm vsto_f ../cli/levo_atm config_vsto_f.xml
arm vsto_d ../cli/stoo_atm config_vsto_d.xml ATM_RH_STORM_LAT=48
wait_arms
for t in a c b f d; do sed -i 's#output_vsto_[a-g]/#OUT/#' output_vsto_$t/RUN_CONFIG.txt; done
cmp_dirs vsto_a vsto_c "A  new clean == old clean"
cmp_dirs vsto_b vsto_f "B  new WB == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  RH_STORM_(LAT|WIDTH)=[^ ]*//g" output_vsto_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vsto_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vsto_d vsto_b "C  CONTROL + RH_STORM_LAT=48 vs 55 -- MUST differ"
touch VSTO_DONE
