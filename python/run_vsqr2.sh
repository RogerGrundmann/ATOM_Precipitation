#!/bin/bash
# 2026-09-29: re-run of run_vsqr0928.sh (killed 09-28 at the user request; its partial output_vsqr_* dirs are kept).
# 2026-09-28 evening: ATM_SEAM_Q_CONSERVE default back to 0 -- -O0 byte check both directions, 1 thread, nm 20 from scratch.
# old = cli/sqf_new_atm (32b949c -O0, default 2), new = HEAD -O0.
#   A  new clean == old SEAM_Q_CONSERVE=0;  B  new SEAM_Q_CONSERVE=2 == old clean;  C  new clean != old clean
set -u; cd "$(dirname "$0")"; rm -f VSQR2_DONE
. ./verify_lib.sh
[ -e ../cli/sqr_new_atm ] || ./build_o0.sh 113221f sqr_new || { touch VSQR2_DONE; exit 1; }
for t in a b c d; do mkdir output_vsqr2_$t || { touch VSQR2_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vsqr2_$t/#" config_vconv_new.xml > config_vsqr2_$t.xml; done
arm vsqr2_a ../cli/sqr_new_atm config_vsqr2_a.xml
arm vsqr2_b ../cli/sqf_new_atm config_vsqr2_b.xml ATM_SEAM_Q_CONSERVE=0
arm vsqr2_c ../cli/sqr_new_atm config_vsqr2_c.xml ATM_SEAM_Q_CONSERVE=2
arm vsqr2_d ../cli/sqf_new_atm config_vsqr2_d.xml
wait_arms
for t in a b c d; do sed -i 's#output_vsqr2_[a-d]/#OUT/#' output_vsqr2_$t/RUN_CONFIG.txt; done
cmp_dirs vsqr2_a vsqr2_b "A  FLIP == old+knob"
diff output_vsqr2_a/RUN_CONFIG.txt output_vsqr2_b/RUN_CONFIG.txt | grep -o 'SEAM_Q_CONSERVE=[^ ]*'
cmp_dirs vsqr2_c vsqr2_d "B  new+knob=2 == old clean"
diff output_vsqr2_c/RUN_CONFIG.txt output_vsqr2_d/RUN_CONFIG.txt | grep -o 'SEAM_Q_CONSERVE=[^ ]*'
want_differ vsqr2_a vsqr2_d "C  CONTROL new clean vs old clean -- MUST differ"
touch VSQR2_DONE
