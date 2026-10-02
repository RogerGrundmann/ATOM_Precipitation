#!/bin/bash
# 2026-10-02: ATM_RH_LAND=2 -- -O0 byte check, 1 thread, nm 20 from scratch. old = cli/rgo_atm (8843fde), new = cli/rl2o_atm.
#   A new clean == old clean;  B new RH_LAND=1 == old RH_LAND=1 (mode 1 unchanged);  C CONTROL new RH_LAND=2 vs old RH_LAND=1, RH_OCEAN 0.82 both -- MUST differ
set -u; cd "$(dirname "$0")"; rm -f VRL2_DONE
. ./verify_lib.sh
for t in a b c d e f; do mkdir output_vrl2_$t || { touch VRL2_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vrl2_$t/#" config_vconv_new.xml > config_vrl2_$t.xml; done
arm vrl2_a ../cli/rl2o_atm config_vrl2_a.xml
arm vrl2_b ../cli/rgo_atm  config_vrl2_b.xml
arm vrl2_c ../cli/rl2o_atm config_vrl2_c.xml ATM_RH_LAND=1
arm vrl2_d ../cli/rgo_atm  config_vrl2_d.xml ATM_RH_LAND=1
arm vrl2_e ../cli/rl2o_atm config_vrl2_e.xml ATM_RH_LAND=2 ATM_RH_OCEAN=0.82
arm vrl2_f ../cli/rgo_atm  config_vrl2_f.xml ATM_RH_LAND=1 ATM_RH_OCEAN=0.82
wait_arms
for t in a b c d e f; do sed -i 's#output_vrl2_[a-f]/#OUT/#' output_vrl2_$t/RUN_CONFIG.txt; done
cmp_dirs vrl2_a vrl2_b "A  new clean == old clean"
cmp_dirs vrl2_c vrl2_d "B  new RH_LAND=1 == old RH_LAND=1"
want_differ vrl2_e vrl2_f "C  CONTROL new RH_LAND=2 vs old RH_LAND=1 (RH_OCEAN 0.82) -- MUST differ"
touch VRL2_DONE
