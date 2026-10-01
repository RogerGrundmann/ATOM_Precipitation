#!/bin/bash
# 2026-10-01: re-run of run_vmvc.sh's arm c, which pointed at a non-existent binary (cli/mvc_atm_OLD, a sed slip;
# the script is corrected). old binary cli/mtv_atm clean, 1 thread, nm 20 from scratch, into output_vmvc_c2; then the
# comparisons A, B, C1 against it.
set -u; cd "$(dirname "$0")"; rm -f VMVC_C2_DONE
. ./verify_lib.sh
mkdir output_vmvc_c2 || { touch VMVC_C2_DONE; exit 1; }
sed "s#output_vconv_new/#output_vmvc_c2/#" config_vconv_new.xml > config_vmvc_c2.xml
arm vmvc_c2 ../cli/mtv_atm config_vmvc_c2.xml
wait_arms
while [ ! -e VMVC_DONE ]; do sleep 10; done          # the a/b/d arms of run_vmvc.sh
sed -i 's#output_vmvc_c2/#OUT/#' output_vmvc_c2/RUN_CONFIG.txt
cmp_dirs vmvc_a vmvc_c2 "A  new clean == old clean"
cmp_dirs vmvc_b vmvc_c2 "B  new VEL_COEFF=0 == old clean"
for t in a b; do diff <(sed 's/  MC_VEL_COEFF=0\*\?//' output_vmvc_$t/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vmvc_c2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$t: RUN_CONFIG differs only by the banner token" || echo "$t: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vmvc_d vmvc_c2 "C1 CONTROL new VEL_COEFF=1 vs old clean -- MUST differ"
touch VMVC_C2_DONE
