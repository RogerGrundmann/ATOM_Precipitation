#!/bin/bash
# 2026-10-01: ATM_MC_VEL_COEFF (MC-TV), new knob default 0 -- -O0 byte check, 1 thread, nm 20 from scratch.
# old = cli/mtv_atm (01d5dc8 + ATM_MC_T_COEFF + ATM_MC_UV_DETRAIN), new = cli/mvc_atm (+ ATM_MC_VEL_COEFF).
#   A  new clean == old clean;  B  new VEL_COEFF=0 == old clean;  C1  CONTROL new VEL_COEFF=1 != old;
#   C2  new all three =1 == old with T_COEFF=1 UV_DETRAIN=1 must DIFFER (VEL_COEFF acts on top of the pair)
set -u; cd "$(dirname "$0")"; rm -f VMVC_DONE
. ./verify_lib.sh
for t in a b c d e f; do mkdir output_vmvc_$t || { touch VMVC_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vmvc_$t/#" config_vconv_new.xml > config_vmvc_$t.xml; done
arm vmvc_a ../cli/mvc_atm config_vmvc_a.xml
arm vmvc_b ../cli/mvc_atm config_vmvc_b.xml ATM_MC_VEL_COEFF=0
arm vmvc_c ../cli/mtv_atm config_vmvc_c.xml
arm vmvc_d ../cli/mvc_atm config_vmvc_d.xml ATM_MC_VEL_COEFF=1
arm vmvc_e ../cli/mvc_atm config_vmvc_e.xml ATM_MC_T_COEFF=1 ATM_MC_UV_DETRAIN=1 ATM_MC_VEL_COEFF=1
arm vmvc_f ../cli/mtv_atm config_vmvc_f.xml ATM_MC_T_COEFF=1 ATM_MC_UV_DETRAIN=1
wait_arms
for t in a b c d e f; do sed -i 's#output_vmvc_[a-e]/#OUT/#' output_vmvc_$t/RUN_CONFIG.txt; done
cmp_dirs vmvc_a vmvc_c "A  new clean == old clean"
cmp_dirs vmvc_b vmvc_c "B  new both=0 == old clean"
for t in a b; do diff <(sed 's/  MC_VEL_COEFF=0\*\?//' output_vmvc_$t/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vmvc_c/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$t: RUN_CONFIG differs only by the banner token" || echo "$t: RUN_CONFIG differs BEYOND the banner tokens"; done
want_differ vmvc_d vmvc_c "C1 CONTROL new VEL_COEFF=1 vs old clean -- MUST differ"
want_differ vmvc_e vmvc_f "C2 CONTROL new all three vs old T_COEFF+UV_DETRAIN -- MUST differ"
touch VMVC_DONE
