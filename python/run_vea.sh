#!/bin/bash
# 2026-10-03: ATM_RH_LAND_EAST (default 0) -- -O0 byte check, 1 thread, nm 20 from scratch. Waits for the vss check.
# old = cli/sso_atm (9aba340 + ATM_RH_SIGMA_SFC), new = cli/eao_atm (+ ATM_RH_LAND_EAST). The knob lives inside the ATM_RH_LAND path, which the
# default never runs, so the off branch is checked twice: on the default (A) and with ATM_RH_LAND=2 (B).
#   A  new clean == old clean;  B  new RH_LAND=2 EAST=0 == old RH_LAND=2;  C  CONTROL new RH_LAND=2 EAST=1 != old RH_LAND=2
set -u; cd "$(dirname "$0")"; rm -f VEA_DONE
until [ -e VSS_DONE ]; do sleep 20; done
. ./verify_lib.sh
for t in a b c d f; do mkdir output_vea_$t || { touch VEA_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vea_$t/#" config_vconv_new.xml > config_vea_$t.xml; done
arm vea_a ../cli/eao_atm config_vea_a.xml
arm vea_c ../cli/sso_atm config_vea_c.xml
arm vea_b ../cli/eao_atm config_vea_b.xml ATM_RH_LAND=2 ATM_RH_LAND_EAST=0
arm vea_f ../cli/sso_atm config_vea_f.xml ATM_RH_LAND=2
arm vea_d ../cli/eao_atm config_vea_d.xml ATM_RH_LAND=2 ATM_RH_LAND_EAST=1
wait_arms
for t in a b c d f; do sed -i 's#output_vea_[a-f]/#OUT/#' output_vea_$t/RUN_CONFIG.txt; done
cmp_dirs vea_a vea_c "A  new clean == old clean"
cmp_dirs vea_b vea_f "B  new RH_LAND=2 EAST=0 == old RH_LAND=2"
for p in "a c" "b f"; do set -- $p; diff <(sed -E "s/  RH_LAND_EAST=[^ ]*//g" output_vea_$1/RUN_CONFIG.txt | tr '\n' ' ' | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') \
                  <(tr '\n' ' ' < output_vea_$2/RUN_CONFIG.txt | sed 's/ *AGCM: \[RUN CONFIG\] knobs: */ /g' | tr -s ' ') >/dev/null \
  && echo "$1: RUN_CONFIG differs only by the banner token" || echo "$1: RUN_CONFIG differs BEYOND the banner token"; done
want_differ vea_d vea_f "C  CONTROL new RH_LAND=2 EAST=1 vs old RH_LAND=2 -- MUST differ"
grep -a "RH-LAND\]\|RH-LAND-EAST" vea_d.log vea_f.log
touch VEA_DONE
