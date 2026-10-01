#!/bin/bash
# 2026-10-01: ATM_MC_T_COEFF + ATM_MC_UV_DETRAIN (MC-TV), two new knobs default 0 -- -O0 byte check, 1 thread, nm 20
# from scratch. old = cli/qdf_atm (3f82fe5 + ATM_Q_DIFF_FLUX), new = cli/mtv_atm (+ print-only f9c6b5f/01d5dc8, off here,
# + the two knobs).
#   A  new clean == old clean;  B  new both =0 == old clean;  C1/C2  CONTROLS new T_COEFF=1 / UV_DETRAIN=1 != old clean
set -u; cd "$(dirname "$0")"; rm -f VMTV_DONE
. ./verify_lib.sh
for t in a b c d e; do mkdir output_vmtv_$t || { touch VMTV_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vmtv_$t/#" config_vconv_new.xml > config_vmtv_$t.xml; done
arm vmtv_a ../cli/mtv_atm config_vmtv_a.xml
arm vmtv_b ../cli/mtv_atm config_vmtv_b.xml ATM_MC_T_COEFF=0 ATM_MC_UV_DETRAIN=0
arm vmtv_c ../cli/qdf_atm config_vmtv_c.xml
arm vmtv_d ../cli/mtv_atm config_vmtv_d.xml ATM_MC_T_COEFF=1
arm vmtv_e ../cli/mtv_atm config_vmtv_e.xml ATM_MC_UV_DETRAIN=1
wait_arms
for t in a b c d e; do sed -i 's#output_vmtv_[a-e]/#OUT/#' output_vmtv_$t/RUN_CONFIG.txt; done
cmp_dirs vmtv_a vmtv_c "A  new clean == old clean"
cmp_dirs vmtv_b vmtv_c "B  new both=0 == old clean"
for t in a b; do diff <(sed 's/  MC_T_COEFF=0\*\?//; s/  MC_UV_DETRAIN=0\*\?//' output_vmtv_$t/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vmtv_c/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$t: RUN_CONFIG differs only by the two banner tokens" || echo "$t: RUN_CONFIG differs BEYOND the banner tokens"; done
want_differ vmtv_d vmtv_c "C1 CONTROL new T_COEFF=1 vs old clean -- MUST differ"
want_differ vmtv_e vmtv_c "C2 CONTROL new UV_DETRAIN=1 vs old clean -- MUST differ"
touch VMTV_DONE
