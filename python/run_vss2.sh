#!/bin/bash
# 2026-10-06: ATM_RH_STORM_SST (+ _REF, _MAX; default 0 = off) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/stoo_atm (4fb4b40 source, -O0), new = cli/ssso_atm (HEAD + working tree, -O0). WB = python/working_branch.env (wb59 stack).
#   A  new clean == old clean      B  new WB == old WB      C  CONTROL new WB + ATM_RH_STORM_SST=0.007 != new WB
# SUPERSEDED by run_vss3.sh (ATM_RH_STORM_SST_LAT was added after this check was launched).
set -u; cd "$(dirname "$0")"; rm -f VSS2_DONE
../python/build_o0.sh WORKTREE ssso || { touch VSS2_DONE; exit 1; }
. ./verify_lib.sh
for t in a c b f d; do mkdir output_vss2_$t || { touch VSS2_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vss2_$t/#" config_vconv_new.xml > config_vss2_$t.xml; done
arm vss2_a ../cli/ssso_atm config_vss2_a.xml
arm vss2_c ../cli/stoo_atm config_vss2_c.xml
. ./working_branch.env
arm vss2_b ../cli/ssso_atm config_vss2_b.xml
arm vss2_f ../cli/stoo_atm config_vss2_f.xml
arm vss2_d ../cli/ssso_atm config_vss2_d.xml ATM_RH_STORM_SST=0.007
wait_arms
for t in a c b f d; do sed -i 's#output_vss2_[a-g]/#OUT/#' output_vss2_$t/RUN_CONFIG.txt; done
cmp_dirs vss2_a vss2_c "A  new clean == old clean"
cmp_dirs vss2_b vss2_f "B  new WB == old WB"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  RH_STORM_SST(_REF|_MAX)?=[^ ]*//g" output_vss2_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vss2_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vss2_d vss2_b "C  CONTROL + RH_STORM_SST=0.007 vs off -- MUST differ"
touch VSS2_DONE
