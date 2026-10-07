#!/bin/bash
# 2026-10-07: ATM_RH_OCEAN_ML_LAT0 / ATM_RH_OCEAN_ML_LAT1 (defaults 30 / 40 = shipped) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/ss3o_atm (2e330a2 source = HEAD source, -O0), new = cli/mlto_atm (HEAD + the two knobs, -O0). WB = python/working_branch.env (wb66 stack,
# ATM_RH_OCEAN_ML=1500, so arm B runs the rewritten taper expression).
#   A  new clean == old clean      B  new WB == old WB      C  CONTROL new WB + ATM_RH_OCEAN_ML_LAT0=22 != new WB
# RESULT (2026-10-07 09:22): A and B identical 13 of 14 (RUN_CONFIG.txt, banner tokens only); control C differs 13 of 14 -> PASS.
set -u; cd "$(dirname "$0")"; rm -f VMLT_DONE
../python/build_o0.sh WORKTREE mlto || { touch VMLT_DONE; exit 1; }
. ./verify_lib.sh
for t in a c b f d; do mkdir output_vmlt_$t || { touch VMLT_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vmlt_$t/#" config_vconv_new.xml > config_vmlt_$t.xml; done
arm vmlt_a ../cli/mlto_atm config_vmlt_a.xml
arm vmlt_c ../cli/ss3o_atm config_vmlt_c.xml
. ./working_branch.env
arm vmlt_b ../cli/mlto_atm config_vmlt_b.xml
arm vmlt_f ../cli/ss3o_atm config_vmlt_f.xml
arm vmlt_d ../cli/mlto_atm config_vmlt_d.xml ATM_RH_OCEAN_ML_LAT0=22
wait_arms
for t in a c b f d; do sed -i 's#output_vmlt_[a-g]/#OUT/#' output_vmlt_$t/RUN_CONFIG.txt; done
cmp_dirs vmlt_a vmlt_c "A  new clean == old clean"
cmp_dirs vmlt_b vmlt_f "B  new WB == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  RH_OCEAN_ML_LAT[01]=[^ ]*//g" output_vmlt_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vmlt_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vmlt_d vmlt_b "C  CONTROL + RH_OCEAN_ML_LAT0=22 vs 30 -- MUST differ"
touch VMLT_DONE
