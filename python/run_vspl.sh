#!/bin/bash
# 2026-10-06: ATM_RH_STORM_POLAR (default 1.0 = off) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/evwo_atm (7c9c2d7 source, -O0; no source change since), new = cli/splo_atm (HEAD + working tree, -O0). WB = python/working_branch.env (wb49 stack).
#   A  new clean == old clean      B  new WB == old WB      C  CONTROL new WB + ATM_RH_STORM_POLAR=1.08 != new WB
# RESULT (11:14): A and B identical 13 of 14 (RUN_CONFIG.txt, banner token only); control C differs 12 of 14 -> PASS.
set -u; cd "$(dirname "$0")"; rm -f VSPL_DONE
../python/build_o0.sh WORKTREE splo || { touch VSPL_DONE; exit 1; }
. ./verify_lib.sh
for t in a c b f d; do mkdir output_vspl_$t || { touch VSPL_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vspl_$t/#" config_vconv_new.xml > config_vspl_$t.xml; done
arm vspl_a ../cli/splo_atm config_vspl_a.xml
arm vspl_c ../cli/evwo_atm config_vspl_c.xml
. ./working_branch.env
arm vspl_b ../cli/splo_atm config_vspl_b.xml
arm vspl_f ../cli/evwo_atm config_vspl_f.xml
arm vspl_d ../cli/splo_atm config_vspl_d.xml ATM_RH_STORM_POLAR=1.08
wait_arms
for t in a c b f d; do sed -i 's#output_vspl_[a-g]/#OUT/#' output_vspl_$t/RUN_CONFIG.txt; done
cmp_dirs vspl_a vspl_c "A  new clean == old clean"
cmp_dirs vspl_b vspl_f "B  new WB == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  RH_STORM_POLAR=[^ ]*//g" output_vspl_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vspl_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vspl_d vspl_b "C  CONTROL + RH_STORM_POLAR=1.08 vs off -- MUST differ"
touch VSPL_DONE
