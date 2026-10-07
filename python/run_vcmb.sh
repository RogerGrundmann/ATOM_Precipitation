#!/bin/bash
# 2026-10-07: ATM_MC_CMB_OCEAN (default 1 = off) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/mlto_atm (9d2e3ec source = HEAD source, -O0), new = cli/cmbo_atm (HEAD + the knob, -O0). WB = python/working_branch.env (wb68 stack).
#   A  new clean == old clean      B  new WB == old WB      C  CONTROL new WB + ATM_MC_CMB_OCEAN=0.65 != new WB
# RESULT (2026-10-07 10:23): A and B identical 13 of 14 (RUN_CONFIG.txt, banner token only); control C differs 13 of 14 -> PASS.
set -u; cd "$(dirname "$0")"; rm -f VCMB_DONE
../python/build_o0.sh WORKTREE cmbo || { touch VCMB_DONE; exit 1; }
. ./verify_lib.sh
for t in a c b f d; do mkdir output_vcmb_$t || { touch VCMB_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vcmb_$t/#" config_vconv_new.xml > config_vcmb_$t.xml; done
arm vcmb_a ../cli/cmbo_atm config_vcmb_a.xml
arm vcmb_c ../cli/mlto_atm config_vcmb_c.xml
. ./working_branch.env
arm vcmb_b ../cli/cmbo_atm config_vcmb_b.xml
arm vcmb_f ../cli/mlto_atm config_vcmb_f.xml
arm vcmb_d ../cli/cmbo_atm config_vcmb_d.xml ATM_MC_CMB_OCEAN=0.65
wait_arms
for t in a c b f d; do sed -i 's#output_vcmb_[a-g]/#OUT/#' output_vcmb_$t/RUN_CONFIG.txt; done
cmp_dirs vcmb_a vcmb_c "A  new clean == old clean"
cmp_dirs vcmb_b vcmb_f "B  new WB == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  MC_CMB_OCEAN=[^ ]*//g" output_vcmb_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vcmb_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vcmb_d vcmb_b "C  CONTROL + MC_CMB_OCEAN=0.65 vs 1 -- MUST differ"
touch VCMB_DONE
