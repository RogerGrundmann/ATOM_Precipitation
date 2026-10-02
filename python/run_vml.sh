#!/bin/bash
# 2026-10-02: ATM_MC_T_ADD, ATM_MC_Q_ADD (=shipped constants), ATM_MC_ML_PARCEL (default 0) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/pbo_atm (5b6c8ba), new = cli/mlo_atm (+ the knob).
#   A  new clean == old clean;  B  new MC_knobs=0 == old clean;  C  CONTROL new MC_knobs=1 != old clean
set -u; cd "$(dirname "$0")"; rm -f VML_DONE
. ./verify_lib.sh
for t in a b c d; do mkdir output_vml_$t || { touch VML_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vml_$t/#" config_vconv_new.xml > config_vml_$t.xml; done
arm vml_a ../cli/mlo_atm config_vml_a.xml
arm vml_b ../cli/mlo_atm config_vml_b.xml ATM_MC_ML_PARCEL=0 ATM_MC_T_ADD=0.2 ATM_MC_Q_ADD=1.0e-4
arm vml_c ../cli/pbo_atm  config_vml_c.xml
arm vml_d ../cli/mlo_atm config_vml_d.xml ATM_MC_T_ADD=0.4
wait_arms
for t in a b c d; do sed -i 's#output_vml_[a-d]/#OUT/#' output_vml_$t/RUN_CONFIG.txt; done
cmp_dirs vml_a vml_c "A  new clean == old clean"
cmp_dirs vml_b vml_c "B  new MC_knobs=0 == old clean"
for t in a b; do diff <(sed -E "s/  MC_(ML_PARCEL|Q_ADD|T_ADD)=[^ ]*//g" output_vml_$t/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vml_c/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$t: RUN_CONFIG differs only by the banner token" || echo "$t: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vml_d vml_c "C  CONTROL new MC_knobs=1 vs old clean -- MUST differ"
touch VML_DONE
