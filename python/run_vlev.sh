#!/bin/bash
# 2026-10-06: ATM_LAND_EVAP (default 0 = off) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/splo_atm (38d5095 source, -O0; no source change since), new = cli/levo_atm (HEAD + working tree, -O0). WB = python/working_branch.env (wb57 stack,
# ATM_WATER_CLOSURE=1 and the evaporation knobs on, so arm B runs the flux branch).
#   A  new clean == old clean      B  new WB == old WB      C  CONTROL new WB + ATM_LAND_EVAP=1 != new WB
# RESULT (12:07): A and B identical 13 of 14 (RUN_CONFIG.txt, banner token only); control C differs 12 of 14 -> PASS.
set -u; cd "$(dirname "$0")"; rm -f VLEV_DONE
../python/build_o0.sh WORKTREE levo || { touch VLEV_DONE; exit 1; }
. ./verify_lib.sh
for t in a c b f d; do mkdir output_vlev_$t || { touch VLEV_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vlev_$t/#" config_vconv_new.xml > config_vlev_$t.xml; done
arm vlev_a ../cli/levo_atm config_vlev_a.xml
arm vlev_c ../cli/splo_atm config_vlev_c.xml
. ./working_branch.env
arm vlev_b ../cli/levo_atm config_vlev_b.xml
arm vlev_f ../cli/splo_atm config_vlev_f.xml
arm vlev_d ../cli/levo_atm config_vlev_d.xml ATM_LAND_EVAP=1
wait_arms
for t in a c b f d; do sed -i 's#output_vlev_[a-g]/#OUT/#' output_vlev_$t/RUN_CONFIG.txt; done
cmp_dirs vlev_a vlev_c "A  new clean == old clean"
cmp_dirs vlev_b vlev_f "B  new WB == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  LAND_EVAP=[^ ]*//g" output_vlev_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vlev_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vlev_d vlev_b "C  CONTROL + LAND_EVAP=1 vs off -- MUST differ"
touch VLEV_DONE
