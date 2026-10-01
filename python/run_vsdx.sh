#!/bin/bash
# 2026-10-01: ATM_SNOW_DEP_FLUX (SNOW-SUBL), new knob default 0 -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/sdo_atm (7d90e72, no knob), new = cli/sdn_atm (+ knob).
#   A  new clean == old clean;  B  new SNOW_DEP_FLUX=0 == old clean;  C  CONTROL new =1 != old clean
set -u; cd "$(dirname "$0")"; rm -f VSDX_DONE
. ./verify_lib.sh
for t in a b c d; do mkdir output_vsdx_$t || { touch VSDX_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vsdx_$t/#" config_vconv_new.xml > config_vsdx_$t.xml; done
arm vsdx_a ../cli/sdn_atm config_vsdx_a.xml
arm vsdx_b ../cli/sdn_atm config_vsdx_b.xml ATM_SNOW_DEP_FLUX=0
arm vsdx_c ../cli/sdo_atm config_vsdx_c.xml
arm vsdx_d ../cli/sdn_atm config_vsdx_d.xml ATM_SNOW_DEP_FLUX=1
wait_arms
for t in a b c d; do sed -i 's#output_vsdx_[a-d]/#OUT/#' output_vsdx_$t/RUN_CONFIG.txt; done
cmp_dirs vsdx_a vsdx_c "A  new clean == old clean"
diff output_vsdx_a/RUN_CONFIG.txt output_vsdx_c/RUN_CONFIG.txt | grep -o 'SNOW_DEP_FLUX=[^ ]*'
cmp_dirs vsdx_b vsdx_c "B  new SNOW_DEP_FLUX=0 == old clean"
diff output_vsdx_b/RUN_CONFIG.txt output_vsdx_c/RUN_CONFIG.txt | grep -o 'SNOW_DEP_FLUX=[^ ]*'
want_differ vsdx_d vsdx_c "C  CONTROL new =1 vs old clean -- MUST differ"
touch VSDX_DONE
