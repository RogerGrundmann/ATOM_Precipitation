#!/bin/bash
# ATM_SATADJ_FADE both-directions byte check. 1 thread, nm = 20 FROM SCRATCH, arms in parallel.
# From scratch because load_state() would overwrite the condensate the fade acts on.
#   vf_new  new binary, knob unset      -- must equal vf_old
#   vf_old  cli/atm_mn, pre-knob binary
#   vf_on   new binary, ATM_SATADJ_FADE=1  -- the want_differ CONTROL
#   vf_no   new binary, ATM_SATADJ_FADE=2  -- the second mode, must also differ
set -u; cd "$(dirname "$0")"; . ./verify_lib.sh
rm -rf output_vf_new output_vf_old output_vf_on output_vf_no; mkdir -p output_vf_new output_vf_old output_vf_on output_vf_no
arm vf_new ../cli/atm    config_vf_new.xml
arm vf_old ../cli/atm_mn config_vf_old.xml
arm vf_on  ../cli/atm    config_vf_on.xml  ATM_SATADJ_FADE=1
arm vf_no  ../cli/atm    config_vf_no.xml  ATM_SATADJ_FADE=2
wait_arms
cmp_dirs    vf_new vf_old "A  off-branch (new unset vs pre-knob binary)"
want_differ vf_on  vf_new "C1 CONTROL mode 1 (freeze) vs off -- MUST differ"
want_differ vf_no  vf_new "C2 CONTROL mode 2 (no fade) vs off -- MUST differ"
touch VFADE_DONE
