#!/bin/bash
# 2026-09-28: ATM_SEAM_Q_CONSERVE default 2 (32b949c) -- -O0 byte check both directions, 1 thread, nm 20 from scratch.
# old = cli/sq2_new_atm (caf17b7 + mode 2 = a50f497 atm code, but WITHOUT the OROG flip cc4ab02 -> old arms set ATM_OROG_Q_MASS=1), new = 32b949c -O0.
#   A  new clean == old ATM_SEAM_Q_CONSERVE=2;  B  new ATM_SEAM_Q_CONSERVE=0 == old clean;  C  new clean != old clean
set -u; cd "$(dirname "$0")"; rm -f VSQF0928_DONE
. ./verify_lib.sh
[ -e ../cli/sqf_new_atm ] || ./build_o0.sh 32b949c sqf_new || { touch VSQF0928_DONE; exit 1; }
for t in a b c d; do mkdir output_vsqf_$t || { touch VSQF0928_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vsqf_$t/#" config_vconv_new.xml > config_vsqf_$t.xml; done
arm vsqf_a ../cli/sqf_new_atm config_vsqf_a.xml
arm vsqf_b ../cli/sq2_new_atm config_vsqf_b.xml ATM_SEAM_Q_CONSERVE=2 ATM_OROG_Q_MASS=1
arm vsqf_c ../cli/sqf_new_atm config_vsqf_c.xml ATM_SEAM_Q_CONSERVE=0
arm vsqf_d ../cli/sq2_new_atm config_vsqf_d.xml ATM_OROG_Q_MASS=1
wait_arms
for t in a b c d; do sed -i 's#output_vsqf_[a-d]/#OUT/#' output_vsqf_$t/RUN_CONFIG.txt; done
cmp_dirs vsqf_a vsqf_b "A  FLIP == old+knob"
diff output_vsqf_a/RUN_CONFIG.txt output_vsqf_b/RUN_CONFIG.txt | grep -o 'SEAM_Q_CONSERVE=[^ ]*'
cmp_dirs vsqf_c vsqf_d "B  new+knob=0 == old clean"
diff output_vsqf_c/RUN_CONFIG.txt output_vsqf_d/RUN_CONFIG.txt | grep -o 'SEAM_Q_CONSERVE=[^ ]*'
want_differ vsqf_a vsqf_d "C  CONTROL new clean vs old clean -- MUST differ"
touch VSQF0928_DONE
