#!/bin/bash
# 2026-10-01: ATM_MC_Q_NDIM (MC-Q-LEAK remainder), new knob default 0 -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/cdn_atm (30455b6 + SNOW_DEP_FLUX + MC_COND_DEBIT), new = cli/qnd_atm (+ ATM_MC_Q_NDIM).
#   A  new clean == old clean;  B  new MC_Q_NDIM=0 == old clean;  C  CONTROL new =1 != old clean
set -u; cd "$(dirname "$0")"; rm -f VQND_DONE
. ./verify_lib.sh
for t in a b c d; do mkdir output_vqnd_$t || { touch VQND_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vqnd_$t/#" config_vconv_new.xml > config_vqnd_$t.xml; done
arm vqnd_a ../cli/qnd_atm config_vqnd_a.xml
arm vqnd_b ../cli/qnd_atm config_vqnd_b.xml ATM_MC_Q_NDIM=0
arm vqnd_c ../cli/cdn_atm config_vqnd_c.xml
arm vqnd_d ../cli/qnd_atm config_vqnd_d.xml ATM_MC_Q_NDIM=1
wait_arms
for t in a b c d; do sed -i 's#output_vqnd_[a-d]/#OUT/#' output_vqnd_$t/RUN_CONFIG.txt; done
cmp_dirs vqnd_a vqnd_c "A  new clean == old clean"
diff output_vqnd_a/RUN_CONFIG.txt output_vqnd_c/RUN_CONFIG.txt | grep -o 'MC_Q_NDIM=[^ ]*'
cmp_dirs vqnd_b vqnd_c "B  new MC_Q_NDIM=0 == old clean"
diff output_vqnd_b/RUN_CONFIG.txt output_vqnd_c/RUN_CONFIG.txt | grep -o 'MC_Q_NDIM=[^ ]*'
want_differ vqnd_d vqnd_c "C  CONTROL new =1 vs old clean -- MUST differ"
touch VQND_DONE
