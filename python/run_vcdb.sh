#!/bin/bash
# 2026-10-01: ATM_MC_COND_DEBIT (MC-Q-LEAK), new knob default 0 -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/sdn_atm (7d90e72 + ATM_SNOW_DEP_FLUX), new = cli/cdn_atm (+ the print-only MC-Q probe/census ebb1eed,
# 30455b6, both behind ATM_MC_DIAG which is off here, + ATM_MC_COND_DEBIT).
#   A  new clean == old clean;  B  new MC_COND_DEBIT=0 == old clean;  C  CONTROL new =1 != old clean
set -u; cd "$(dirname "$0")"; rm -f VCDB_DONE
. ./verify_lib.sh
for t in a b c d; do mkdir output_vcdb_$t || { touch VCDB_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vcdb_$t/#" config_vconv_new.xml > config_vcdb_$t.xml; done
arm vcdb_a ../cli/cdn_atm config_vcdb_a.xml
arm vcdb_b ../cli/cdn_atm config_vcdb_b.xml ATM_MC_COND_DEBIT=0
arm vcdb_c ../cli/sdn_atm config_vcdb_c.xml
arm vcdb_d ../cli/cdn_atm config_vcdb_d.xml ATM_MC_COND_DEBIT=1
wait_arms
for t in a b c d; do sed -i 's#output_vcdb_[a-d]/#OUT/#' output_vcdb_$t/RUN_CONFIG.txt; done
cmp_dirs vcdb_a vcdb_c "A  new clean == old clean"
diff output_vcdb_a/RUN_CONFIG.txt output_vcdb_c/RUN_CONFIG.txt | grep -o 'MC_COND_DEBIT=[^ ]*'
cmp_dirs vcdb_b vcdb_c "B  new MC_COND_DEBIT=0 == old clean"
diff output_vcdb_b/RUN_CONFIG.txt output_vcdb_c/RUN_CONFIG.txt | grep -o 'MC_COND_DEBIT=[^ ]*'
want_differ vcdb_d vcdb_c "C  CONTROL new =1 vs old clean -- MUST differ"
touch VCDB_DONE
