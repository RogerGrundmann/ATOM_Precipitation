#!/bin/bash
# 2026-09-29: ATM_SURF_DRAG_CONSISTENT default 0.0 -> 1.0 (user) -- -O0 byte check both directions, 1 thread, nm 20
# from scratch. old = cli/bsc_atm (19f9a0d's atm source, -O0, default 0.0), new = cli/sdf_atm (+ the flip, -O0).
#   A  new clean == old SURF_DRAG=1.0;  B  new SURF_DRAG=0 == old clean;  C  new clean != old clean
set -u; cd "$(dirname "$0")"; rm -f VSDF_DONE
. ./verify_lib.sh
for t in a b c d; do mkdir output_vsdf_$t || { touch VSDF_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vsdf_$t/#" config_vconv_new.xml > config_vsdf_$t.xml; done
arm vsdf_a ../cli/sdf_atm config_vsdf_a.xml
arm vsdf_b ../cli/bsc_atm config_vsdf_b.xml ATM_SURF_DRAG_CONSISTENT=1.0
arm vsdf_c ../cli/sdf_atm config_vsdf_c.xml ATM_SURF_DRAG_CONSISTENT=0
arm vsdf_d ../cli/bsc_atm config_vsdf_d.xml
wait_arms
for t in a b c d; do sed -i 's#output_vsdf_[a-d]/#OUT/#' output_vsdf_$t/RUN_CONFIG.txt; done
cmp_dirs vsdf_a vsdf_b "A  FLIP == old+knob"
diff output_vsdf_a/RUN_CONFIG.txt output_vsdf_b/RUN_CONFIG.txt | grep -o 'SURF_DRAG_CONSISTENT=[^ ]*'
cmp_dirs vsdf_c vsdf_d "B  new+knob=0 == old clean"
diff output_vsdf_c/RUN_CONFIG.txt output_vsdf_d/RUN_CONFIG.txt | grep -o 'SURF_DRAG_CONSISTENT=[^ ]*'
want_differ vsdf_a vsdf_d "C  CONTROL new clean vs old clean -- MUST differ"
touch VSDF_DONE
