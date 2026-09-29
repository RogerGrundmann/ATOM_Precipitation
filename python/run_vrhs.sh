#!/bin/bash
# 2026-09-29: ATM_RH_STORM (new, default 1.0) -- -O0 off-branch byte check, 1 thread, nm 20 from scratch.
# old = cli/ced_atm (c8f1489 -O0), new = cli/rhs_atm (+ the knob, -O0).
set -u; cd "$(dirname "$0")"; rm -f VRHS_DONE
. ./verify_lib.sh
for t in a b c; do mkdir output_vrhs_$t || { touch VRHS_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vrhs_$t/#" config_vconv_new.xml > config_vrhs_$t.xml; done
arm vrhs_a ../cli/rhs_atm config_vrhs_a.xml
arm vrhs_b ../cli/ced_atm config_vrhs_b.xml
arm vrhs_c ../cli/rhs_atm config_vrhs_c.xml ATM_RH_STORM=1.25
wait_arms
for t in a b c; do sed -i 's#output_vrhs_[a-c]/#OUT/#' output_vrhs_$t/RUN_CONFIG.txt; done
cmp_dirs vrhs_a vrhs_b "A  new clean == old clean"
want_differ vrhs_c vrhs_b "C  CONTROL new RH_STORM=1.25 vs old clean -- MUST differ"
touch VRHS_DONE
