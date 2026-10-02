#!/bin/bash
# 2026-10-02: ATM_MC_DEPTH_RAMP (default 0) -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/wbo_atm (14dca2a), new = cli/dro_atm (+ the knob).
#   A  new clean == old clean;  B  new MC_knobs=0 == old clean;  C  CONTROL new MC_knobs=1 != old clean
set -u; cd "$(dirname "$0")"; rm -f VDR_DONE
. ./verify_lib.sh
for t in a b c d; do mkdir output_vdr_$t || { touch VDR_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vdr_$t/#" config_vconv_new.xml > config_vdr_$t.xml; done
arm vdr_a ../cli/dro_atm config_vdr_a.xml
arm vdr_b ../cli/dro_atm config_vdr_b.xml ATM_MC_DEPTH_RAMP=0
arm vdr_c ../cli/wbo_atm  config_vdr_c.xml
arm vdr_d ../cli/dro_atm config_vdr_d.xml ATM_MC_DEPTH_RAMP=1
wait_arms
for t in a b c d; do sed -i 's#output_vdr_[a-d]/#OUT/#' output_vdr_$t/RUN_CONFIG.txt; done
cmp_dirs vdr_a vdr_c "A  new clean == old clean"
cmp_dirs vdr_b vdr_c "B  new MC_knobs=0 == old clean"
for t in a b; do diff <(sed -E "s/  MC_DEPTH_RAMP=[^ ]*//g" output_vdr_$t/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vdr_c/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$t: RUN_CONFIG differs only by the banner token" || echo "$t: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vdr_d vdr_c "C  CONTROL new MC_knobs=1 vs old clean -- MUST differ"
touch VDR_DONE
