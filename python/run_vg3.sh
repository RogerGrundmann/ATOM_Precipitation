#!/bin/bash
# 2026-10-03: ATM_MC_GATE_BLEND=3 (new mode; default 0) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/fin_atm (957b13c source), new = cli/g3o_atm. K = ATM_MC_ML_PARCEL=1 ATM_MC_GATE_BLEND=1: the existing mode 1 must be unchanged.
#   A  new clean == old clean;  B  new K == old K (mode 1 untouched);  C  CONTROL new ML_PARCEL=1 GATE_BLEND=3 != old K
set -u; cd "$(dirname "$0")"; rm -f VG3_DONE
. ./verify_lib.sh
for t in a b c d f; do mkdir output_vg3_$t || { touch VG3_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vg3_$t/#" config_vconv_new.xml > config_vg3_$t.xml; done
K="ATM_MC_ML_PARCEL=1 ATM_MC_GATE_BLEND=1"
arm vg3_a ../cli/g3o_atm config_vg3_a.xml
arm vg3_c ../cli/fin_atm config_vg3_c.xml
arm vg3_b ../cli/g3o_atm config_vg3_b.xml $K
arm vg3_f ../cli/fin_atm config_vg3_f.xml $K
arm vg3_d ../cli/g3o_atm config_vg3_d.xml ATM_MC_ML_PARCEL=1 ATM_MC_GATE_BLEND=3
wait_arms
for t in a b c d f; do sed -i 's#output_vg3_[a-f]/#OUT/#' output_vg3_$t/RUN_CONFIG.txt; done
cmp_dirs vg3_a vg3_c "A  new clean == old clean"
cmp_dirs vg3_b vg3_f "B  new K == old K (mode 1 untouched)"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  NOTOKEN=[^ ]*//g" output_vg3_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vg3_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vg3_d vg3_f "C  CONTROL new GATE_BLEND=3 vs old GATE_BLEND=1 -- MUST differ"
touch VG3_DONE
