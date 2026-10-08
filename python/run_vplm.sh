#!/bin/bash
# 2026-10-08: ATM_HCRIT_SFC_POLAR (+ _LAT) and ATM_RH_LAND_ML (+ _STRENGTH, _LAT0, _LAT1), all default off -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/cmbo_atm (HEAD 5035541 source, -O0), new = cli/plmo_atm (HEAD + the knobs, -O0). WB = python/working_branch.env (wb68 stack).
#   A  new clean == old clean      B  new WB == old WB      C  CONTROL new WB + ATM_HCRIT_SFC_POLAR=1 != new WB      D  CONTROL new WB + ATM_RH_LAND_ML=1500 != new WB
# RESULT (2026-10-08 08:37): A and B identical 13 of 14 (RUN_CONFIG.txt, banner tokens only); controls C and D differ 12 of 14 -> PASS.
set -u; cd "$(dirname "$0")"; rm -f VPLM_DONE
../python/build_o0.sh WORKTREE plmo || { touch VPLM_DONE; exit 1; }
. ./verify_lib.sh
for t in a c b f d e; do mkdir output_vplm_$t || { touch VPLM_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vplm_$t/#" config_vconv_new.xml > config_vplm_$t.xml; done
arm vplm_a ../cli/plmo_atm config_vplm_a.xml
arm vplm_c ../cli/cmbo_atm config_vplm_c.xml
. ./working_branch.env
arm vplm_b ../cli/plmo_atm config_vplm_b.xml
arm vplm_f ../cli/cmbo_atm config_vplm_f.xml
arm vplm_d ../cli/plmo_atm config_vplm_d.xml ATM_HCRIT_SFC_POLAR=1
arm vplm_e ../cli/plmo_atm config_vplm_e.xml ATM_RH_LAND_ML=1500
wait_arms
for t in a c b f d e; do sed -i 's#output_vplm_[a-g]/#OUT/#' output_vplm_$t/RUN_CONFIG.txt; done
cmp_dirs vplm_a vplm_c "A  new clean == old clean"
cmp_dirs vplm_b vplm_f "B  new WB == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  (HCRIT_SFC_POLAR|HCRIT_SFC_POLAR_LAT|RH_LAND_ML|RH_LAND_ML_STRENGTH|RH_LAND_ML_LAT0|RH_LAND_ML_LAT1)=[^ ]*//g" output_vplm_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vplm_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner tokens" || echo "$1: RUN_CONFIG differs BEYOND the banner tokens"; done
want_differ vplm_d vplm_b "C  CONTROL + HCRIT_SFC_POLAR=1 vs 0 -- MUST differ"
want_differ vplm_e vplm_b "D  CONTROL + RH_LAND_ML=1500 vs 0 -- MUST differ"
touch VPLM_DONE
