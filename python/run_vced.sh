#!/bin/bash
# 2026-09-29: ATM_MC_ED_AREA (new, default 0) -- -O0 off-branch byte check, 1 thread, nm 20 from scratch.
# old = cli/sdf_atm (55489e9 -O0), new = cli/ced_atm (+ the knob, -O0).
#   A  new clean == old clean;  C  new ED_AREA=1 != old clean (the knob must fire)
set -u; cd "$(dirname "$0")"; rm -f VCED_DONE
. ./verify_lib.sh
for t in a b c; do mkdir output_vced_$t || { touch VCED_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vced_$t/#" config_vconv_new.xml > config_vced_$t.xml; done
arm vced_a ../cli/ced_atm config_vced_a.xml
arm vced_b ../cli/sdf_atm config_vced_b.xml
arm vced_c ../cli/ced_atm config_vced_c.xml ATM_MC_ED_AREA=1
wait_arms
for t in a b c; do sed -i 's#output_vced_[a-c]/#OUT/#' output_vced_$t/RUN_CONFIG.txt; done
cmp_dirs vced_a vced_b "A  new clean == old clean"
diff output_vced_a/RUN_CONFIG.txt output_vced_b/RUN_CONFIG.txt | head -4 | cut -c1-80
want_differ vced_c vced_b "C  CONTROL new ED_AREA=1 vs old clean -- MUST differ"
touch VCED_DONE
