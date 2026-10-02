#!/bin/bash
# 2026-10-02: ATM_TEQ_WTG, new knob default 0 -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/lo_atm (d2fa7ed), new = cli/wtg_atm (+ the knob).
#   A  new clean == old clean;  B  new TEQ_WTG=0 == old clean;  C  CONTROL new TEQ_WTG=1 != old clean
set -u; cd "$(dirname "$0")"; rm -f VWTG_DONE
. ./verify_lib.sh
for t in a b c d; do mkdir output_vwtg_$t || { touch VWTG_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vwtg_$t/#" config_vconv_new.xml > config_vwtg_$t.xml; done
arm vwtg_a ../cli/wtg_atm config_vwtg_a.xml
arm vwtg_b ../cli/wtg_atm config_vwtg_b.xml ATM_TEQ_WTG=0
arm vwtg_c ../cli/lo_atm  config_vwtg_c.xml
arm vwtg_d ../cli/wtg_atm config_vwtg_d.xml ATM_TEQ_WTG=1
wait_arms
for t in a b c d; do sed -i 's#output_vwtg_[a-d]/#OUT/#' output_vwtg_$t/RUN_CONFIG.txt; done
cmp_dirs vwtg_a vwtg_c "A  new clean == old clean"
cmp_dirs vwtg_b vwtg_c "B  new TEQ_WTG=0 == old clean"
for t in a b; do diff <(sed 's/  TEQ_WTG=0\(\.0\)\?\*\?//' output_vwtg_$t/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vwtg_c/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$t: RUN_CONFIG differs only by the banner token" || echo "$t: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vwtg_d vwtg_c "C  CONTROL new TEQ_WTG=1 vs old clean -- MUST differ"
grep -a "TEQ-WTG" output_vwtg_d/*.log vwtg_d.log 2>/dev/null | head -2
touch VWTG_DONE
