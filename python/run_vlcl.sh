#!/bin/bash
# 2026-10-02: ATM_MC_ML_LCL (default 0) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/rhoc_atm (3f91cb9), new = cli/lclo_atm (+ the knob).
#   A  new clean == old clean;  B  new MC_knobs=0 == old clean;  C  CONTROL new MC_knobs=1 != old clean
set -u; cd "$(dirname "$0")"; rm -f VLCL_DONE
. ./verify_lib.sh
for t in a b c d; do mkdir output_vlcl_$t || { touch VLCL_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vlcl_$t/#" config_vconv_new.xml > config_vlcl_$t.xml; done
arm vlcl_a ../cli/lclo_atm config_vlcl_a.xml
arm vlcl_b ../cli/lclo_atm config_vlcl_b.xml ATM_MC_ML_LCL=0
arm vlcl_c ../cli/rhoc_atm  config_vlcl_c.xml
arm vlcl_d ../cli/lclo_atm config_vlcl_d.xml ATM_MC_ML_PARCEL=2 ATM_MC_ML_LCL=1
wait_arms
for t in a b c d; do sed -i 's#output_vlcl_[a-d]/#OUT/#' output_vlcl_$t/RUN_CONFIG.txt; done
cmp_dirs vlcl_a vlcl_c "A  new clean == old clean"
cmp_dirs vlcl_b vlcl_c "B  new MC_knobs=0 == old clean"
for t in a b; do diff <(sed -E "s/  MC_ML_LCL=[^ ]*//g" output_vlcl_$t/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vlcl_c/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$t: RUN_CONFIG differs only by the banner token" || echo "$t: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vlcl_d vlcl_c "C  CONTROL new MC_knobs=1 vs old clean -- MUST differ"
touch VLCL_DONE
