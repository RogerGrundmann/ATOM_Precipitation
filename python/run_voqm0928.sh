#!/bin/bash
# 2026-09-28: ATM_OROG_Q_MASS flipped ON (cc4ab02) -- -O0 byte check both directions, 1 thread, nm 20 from scratch.
# old = cli/sq2_old_atm (caf17b7 -O0; b12e7d2/d5d0272 change only the ocean), new = cc4ab02 -O0.
#   A  new clean == old ATM_OROG_Q_MASS=1;  B  new ATM_OROG_Q_MASS=0 == old clean;  C  new clean != old clean
set -u; cd "$(dirname "$0")"; rm -f VOQM0928_DONE
. ./verify_lib.sh
[ -e ../cli/oqm_new_atm ] || ./build_o0.sh cc4ab02 oqm_new || { touch VOQM0928_DONE; exit 1; }
for t in a b c d; do mkdir output_voqm_$t || { touch VOQM0928_DONE; exit 1; }
  sed "s#output_vconv_new/#output_voqm_$t/#" config_vconv_new.xml > config_voqm_$t.xml; done
arm voqm_a ../cli/oqm_new_atm config_voqm_a.xml
arm voqm_b ../cli/sq2_old_atm config_voqm_b.xml ATM_OROG_Q_MASS=1
arm voqm_c ../cli/oqm_new_atm config_voqm_c.xml ATM_OROG_Q_MASS=0
arm voqm_d ../cli/sq2_old_atm config_voqm_d.xml
wait_arms
for t in a b c d; do sed -i 's#output_voqm_[a-d]/#OUT/#' output_voqm_$t/RUN_CONFIG.txt; done
cmp_dirs voqm_a voqm_b "A  FLIP == old+knob"
diff output_voqm_a/RUN_CONFIG.txt output_voqm_b/RUN_CONFIG.txt | grep -o 'OROG_Q_MASS=[^ ]*'
cmp_dirs voqm_c voqm_d "B  new+knob=0 == old clean"
diff output_voqm_c/RUN_CONFIG.txt output_voqm_d/RUN_CONFIG.txt | grep -o 'OROG_Q_MASS=[^ ]*'
want_differ voqm_a voqm_d "C  CONTROL new clean vs old clean -- MUST differ"
touch VOQM0928_DONE
