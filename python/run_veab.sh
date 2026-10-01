#!/bin/bash
# 2026-10-01: ATM_MC_ED_ABOVE_BASE, new knob default 0 -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/qcd_atm (f90c462 source), new = cli/eab_atm (+ the knob).
#   A  new clean == old clean;  B  new ED_ABOVE_BASE=0 == old clean;  C  CONTROL new ED_ABOVE_BASE=1 != old clean
set -u; cd "$(dirname "$0")"; rm -f VEAB_DONE
# NOTE: the control differs in only 2 of 14 files (the 28N slice and the banner): on the shipped default the convective
# rain at 20 iterations is small, so the knob changes little -- it fires, weakly.
. ./verify_lib.sh
for t in a b c d; do mkdir output_veab_$t || { touch VEAB_DONE; exit 1; }
  sed "s#output_vconv_new/#output_veab_$t/#" config_vconv_new.xml > config_veab_$t.xml; done
arm veab_a ../cli/eab_atm config_veab_a.xml
arm veab_b ../cli/eab_atm config_veab_b.xml ATM_MC_ED_ABOVE_BASE=0
arm veab_c ../cli/qcd_atm config_veab_c.xml
arm veab_d ../cli/eab_atm config_veab_d.xml ATM_MC_ED_ABOVE_BASE=1
wait_arms
for t in a b c d; do sed -i 's#output_veab_[a-d]/#OUT/#' output_veab_$t/RUN_CONFIG.txt; done
cmp_dirs veab_a veab_c "A  new clean == old clean"
cmp_dirs veab_b veab_c "B  new ED_ABOVE_BASE=0 == old clean"
for t in a b; do diff <(sed 's/  MC_ED_ABOVE_BASE=0\*\?//' output_veab_$t/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_veab_c/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$t: RUN_CONFIG differs only by the banner token" || echo "$t: RUN_CONFIG differs BEYOND the banner token"; done
want_differ veab_d veab_c "C  CONTROL new ED_ABOVE_BASE=1 vs old clean -- MUST differ"
touch VEAB_DONE
