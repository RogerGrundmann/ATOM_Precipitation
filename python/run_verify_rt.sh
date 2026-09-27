#!/bin/bash
# Byte check for the OneCat/ThreeCat retirement + CLI exit code (2026-09-27). Both -O2, 1 thread, nm = 20 from scratch.
# old = cli/atm_ocfix (f611736), new = cli/atm_ret. Default branch (TwoCat) must be identical.
set -u; cd "$(dirname "$0")"; rm -f RT_VERIFY_DONE
. ./verify_lib.sh
arm rt_old ../cli/atm_ocfix config_rt_old.xml
arm rt_new ../cli/atm_ret   config_rt_new.xml
wait_arms
cmp_dirs rt_new rt_old "A  OFF BRANCH (retirement, default TwoCat)"
touch RT_VERIFY_DONE
