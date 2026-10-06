#!/bin/bash
# 2026-10-06: ATM_RH_STORM_SST (+ _REF, _MAX, _LAT; default 0 = off; supersedes run_vss2.sh, which checked the code before _LAT was added) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/stoo_atm (4fb4b40 source, -O0), new = cli/ss3o_atm (HEAD + working tree, -O0). WB = python/working_branch.env (wb59 stack).
#   A  new clean == old clean      B  new WB == old WB      C  CONTROL new WB + ATM_RH_STORM_SST=0.007 != new WB
# RESULT (about 14:10): A and B identical 13 of 14 (RUN_CONFIG.txt, banner tokens only); control C differs 13 of 14 -> PASS.
set -u; cd "$(dirname "$0")"; rm -f VSS3_DONE
../python/build_o0.sh WORKTREE ss3o || { touch VSS3_DONE; exit 1; }
. ./verify_lib.sh
for t in a c b f d; do mkdir output_vss3_$t || { touch VSS3_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vss3_$t/#" config_vconv_new.xml > config_vss3_$t.xml; done
arm vss3_a ../cli/ss3o_atm config_vss3_a.xml
arm vss3_c ../cli/stoo_atm config_vss3_c.xml
. ./working_branch.env
arm vss3_b ../cli/ss3o_atm config_vss3_b.xml
arm vss3_f ../cli/stoo_atm config_vss3_f.xml
arm vss3_d ../cli/ss3o_atm config_vss3_d.xml ATM_RH_STORM_SST=0.007
wait_arms
for t in a c b f d; do sed -i 's#output_vss3_[a-g]/#OUT/#' output_vss3_$t/RUN_CONFIG.txt; done
cmp_dirs vss3_a vss3_c "A  new clean == old clean"
cmp_dirs vss3_b vss3_f "B  new WB == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  RH_STORM_SST(_REF|_MAX|_LAT)?=[^ ]*//g" output_vss3_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vss3_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vss3_d vss3_b "C  CONTROL + RH_STORM_SST=0.007 vs off -- MUST differ"
touch VSS3_DONE
