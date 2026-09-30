#!/bin/bash
# 2026-09-30: plan D ocean XML arm, re-run without HYD_T_FREEZE_SFC (retired in batch 2 -- the <knobs> guard rightly
# refused it). New = cli/xml2_hyd with 4 metric knobs in <knobs>; old = output_vxml_z (ret_b5_hyd, same values in env).
set -u; cd "$(dirname "$0")"; rm -f VXML3_DONE; . ./verify_lib.sh
arm vxml3_y ../cli/xml2_hyd config_vxml3_y.xml; wait_arms
cmp_dirs vxml3_y vxml_z "D  hyd new with <knobs> == old+env"
grep -a '\[RUN CONFIG\] knobs' vxml3_y.log | grep -o 'METRIC_RADIUS=[^ ]*\|RUN_NEUMANN=[^ ]*' | head -2
touch VXML3_DONE
