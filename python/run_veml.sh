#!/bin/bash
# 2026-10-06: ATM_RH_LAND_EAST_ML (+ _STRENGTH, _T; default 0 = off) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/emxo_atm (8920fd5 source, -O0; no source change since), new = cli/emlo_atm (HEAD + working tree, -O0). WB = python/working_branch.env (wb49 stack).
#   A  new clean == old clean
#   B  new WB == old WB
#   C  CONTROL new WB + ATM_RH_LAND_EAST_ML=1500 ATM_RH_LAND_EAST_ML_STRENGTH=0.6 != new WB   (the mixed layer fires)
# RESULT (10:02): A and B identical 13 of 14 (RUN_CONFIG.txt, banner tokens only); control C differs 13 of 14 -> PASS.
set -u; cd "$(dirname "$0")"; rm -f VEML_DONE
../python/build_o0.sh WORKTREE emlo || { touch VEML_DONE; exit 1; }
. ./verify_lib.sh
for t in a c b f d; do mkdir output_veml_$t || { touch VEML_DONE; exit 1; }
  sed "s#output_vconv_new/#output_veml_$t/#" config_vconv_new.xml > config_veml_$t.xml; done
arm veml_a ../cli/emlo_atm config_veml_a.xml
arm veml_c ../cli/emxo_atm config_veml_c.xml
. ./working_branch.env
arm veml_b ../cli/emlo_atm config_veml_b.xml
arm veml_f ../cli/emxo_atm config_veml_f.xml
arm veml_d ../cli/emlo_atm config_veml_d.xml ATM_RH_LAND_EAST_ML=1500 ATM_RH_LAND_EAST_ML_STRENGTH=0.6
wait_arms
for t in a c b f d; do sed -i 's#output_veml_[a-g]/#OUT/#' output_veml_$t/RUN_CONFIG.txt; done
cmp_dirs veml_a veml_c "A  new clean == old clean"
cmp_dirs veml_b veml_f "B  new WB == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  RH_LAND_EAST_ML(_STRENGTH|_T)?=[^ ]*//g" output_veml_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_veml_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ veml_d veml_b "C  CONTROL + RH_LAND_EAST_ML=1500 x 0.6 vs off -- MUST differ"
touch VEML_DONE
