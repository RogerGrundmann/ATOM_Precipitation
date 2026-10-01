#!/bin/bash
# 2026-10-01: ATM_Q_DIFF_FLUX (DIFF-LEAK), new knob default 0 -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/qnd_atm (8a9e59e + MC_Q_NDIM), new = cli/qdf_atm (+ print-only CWB-TD probes db650eb..3f82fe5, off here, + the knob).
#   A  new clean == old clean;  B  new Q_DIFF_FLUX=0 == old clean;  C  CONTROL new =1 != old clean
set -u; cd "$(dirname "$0")"; rm -f VQDF_DONE
. ./verify_lib.sh
for t in a b c d; do mkdir output_vqdf_$t || { touch VQDF_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vqdf_$t/#" config_vconv_new.xml > config_vqdf_$t.xml; done
arm vqdf_a ../cli/qdf_atm config_vqdf_a.xml
arm vqdf_b ../cli/qdf_atm config_vqdf_b.xml ATM_Q_DIFF_FLUX=0
arm vqdf_c ../cli/qnd_atm config_vqdf_c.xml
arm vqdf_d ../cli/qdf_atm config_vqdf_d.xml ATM_Q_DIFF_FLUX=1
wait_arms
for t in a b c d; do sed -i 's#output_vqdf_[a-d]/#OUT/#' output_vqdf_$t/RUN_CONFIG.txt; done
cmp_dirs vqdf_a vqdf_c "A  new clean == old clean"
diff output_vqdf_a/RUN_CONFIG.txt output_vqdf_c/RUN_CONFIG.txt | grep -o 'Q_DIFF_FLUX=[^ ]*'
cmp_dirs vqdf_b vqdf_c "B  new Q_DIFF_FLUX=0 == old clean"
diff output_vqdf_b/RUN_CONFIG.txt output_vqdf_c/RUN_CONFIG.txt | grep -o 'Q_DIFF_FLUX=[^ ]*'
want_differ vqdf_d vqdf_c "C  CONTROL new =1 vs old clean -- MUST differ"
touch VQDF_DONE
