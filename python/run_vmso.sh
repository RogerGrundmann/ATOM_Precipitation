#!/bin/bash
# 2026-10-04: ATM_MC_MB_SAT_OCEAN (default 0) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/h27o_atm (HEAD 755f191, -O0), new = cli/mso_atm (HEAD + working tree, -O0). WB = python/working_branch.env (the wb27 stack).
#   A  new clean == old clean
#   B  new WB + MB_SAT_OCEAN=0 == old WB
#   C  CONTROL new WB + MB_SAT_OCEAN=0.03 != old WB
#   D  new WB + MB_SAT_LAND=0.04 == old WB + MB_SAT_LAND=0.04   (the land branch of mbSat is unchanged)
set -u; cd "$(dirname "$0")"; rm -f VMSO_DONE
. ./verify_lib.sh
for t in a c b f d e g; do mkdir output_vmso_$t || { touch VMSO_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vmso_$t/#" config_vconv_new.xml > config_vmso_$t.xml; done
arm vmso_a ../cli/mso_atm  config_vmso_a.xml
arm vmso_c ../cli/h27o_atm config_vmso_c.xml
. ./working_branch.env
arm vmso_b ../cli/mso_atm  config_vmso_b.xml ATM_MC_MB_SAT_OCEAN=0
arm vmso_f ../cli/h27o_atm config_vmso_f.xml
arm vmso_d ../cli/mso_atm  config_vmso_d.xml ATM_MC_MB_SAT_OCEAN=0.03
arm vmso_e ../cli/mso_atm  config_vmso_e.xml ATM_MC_MB_SAT_LAND=0.04
arm vmso_g ../cli/h27o_atm config_vmso_g.xml ATM_MC_MB_SAT_LAND=0.04
wait_arms
for t in a c b f d e g; do sed -i 's#output_vmso_[a-g]/#OUT/#' output_vmso_$t/RUN_CONFIG.txt; done
grep -a "nm\b\|nm=" output_vmso_a/RUN_CONFIG.txt | head -2
cmp_dirs vmso_a vmso_c "A  new clean == old clean"
cmp_dirs vmso_b vmso_f "B  new WB MB_SAT_OCEAN=0 == old WB"
cmp_dirs vmso_e vmso_g "D  new WB MB_SAT_LAND=0.04 == old WB MB_SAT_LAND=0.04"
for p in "a c" "b f" "e g"; do set -- $p; diff <(sed -E "s/  MC_MB_SAT_OCEAN=[^ ]*//g" output_vmso_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vmso_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vmso_d vmso_f "C  CONTROL new WB MB_SAT_OCEAN=0.03 vs old WB -- MUST differ"
touch VMSO_DONE
