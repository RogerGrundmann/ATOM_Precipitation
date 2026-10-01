#!/bin/bash
# 2026-10-01: ATM_MC_QC_DETRAIN (MC-TV), new knob default 0 -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/mvc_atm (c56c12b + ATM_MC_VEL_COEFF), new = cli/qcd_atm (+ the print-only [MC-QC] line, off here, + the knob).
#   A  new clean == old clean;  B  new QC_DETRAIN=0 == old clean;  C  CONTROL new QC_DETRAIN=1 != old clean
set -u; cd "$(dirname "$0")"; rm -f VQCD_DONE
. ./verify_lib.sh
for t in a b c d; do mkdir output_vqcd_$t || { touch VQCD_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vqcd_$t/#" config_vconv_new.xml > config_vqcd_$t.xml; done
arm vqcd_a ../cli/qcd_atm config_vqcd_a.xml
arm vqcd_b ../cli/qcd_atm config_vqcd_b.xml ATM_MC_QC_DETRAIN=0
arm vqcd_c ../cli/mvc_atm config_vqcd_c.xml
arm vqcd_d ../cli/qcd_atm config_vqcd_d.xml ATM_MC_QC_DETRAIN=1
wait_arms
for t in a b c d; do sed -i 's#output_vqcd_[a-d]/#OUT/#' output_vqcd_$t/RUN_CONFIG.txt; done
cmp_dirs vqcd_a vqcd_c "A  new clean == old clean"
cmp_dirs vqcd_b vqcd_c "B  new QC_DETRAIN=0 == old clean"
for t in a b; do diff <(sed 's/  MC_QC_DETRAIN=0\*\?//' output_vqcd_$t/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vqcd_c/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$t: RUN_CONFIG differs only by the banner token" || echo "$t: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vqcd_d vqcd_c "C  CONTROL new QC_DETRAIN=1 vs old clean -- MUST differ"
touch VQCD_DONE
