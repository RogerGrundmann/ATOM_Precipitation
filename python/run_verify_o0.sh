#!/bin/bash
# -O0 ATTRIBUTION of the -O2 upwind byte-check failure (2026-09-27). At -O2 the new source (ZeroCat/OneCat/ThreeCat
# dropped + ATM_PRECIP_UPWIND, knob off) differs from the old in 4 S_r cells, all signed zero (-0.00000000 vs
# 0.00000000). Built at -O0 (no floating-point reordering) from HEAD (f611736) and HEAD + the diff: if identical,
# the -O2 difference is compiler reassociation under -ffast-math, not logic.
set -u; cd "$(dirname "$0")"; rm -f O0_VERIFY_DONE
. ./verify_lib.sh
arm o0_old ../cli/atm_o0old config_o0_old.xml
arm o0_new ../cli/atm_o0new config_o0_new.xml
wait_arms
cmp_dirs o0_new o0_old "A  OFF BRANCH at -O0 (new vs old)"
touch O0_VERIFY_DONE
