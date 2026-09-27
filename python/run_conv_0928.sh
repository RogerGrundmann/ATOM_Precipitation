#!/bin/bash
# FOR 2026-09-28 -- LAUNCH BY HAND (NOT QUEUED). Is "the converged rain is ~9 % of NASA" real, and what eats the rain?
# Background (2026-09-27): TwoCat's 3-pass rain iteration does not converge; the model writes pass 3 (0-15 deg: passes
# 370 / 93 / 3465 mm/a). ATM_PRECIP_UPWIND=1 converges exactly (|p3-p2| = 0) and rains 92 mm/a (default branch) /
# 24.5 (closure) at 600 -> 640, against 988 / 831 shipped and NASA 978.
#
# 1) -O0 byte check (both sides -O0, python/build_o0.sh): old = 92f0dfa, new = the working tree/HEAD with
#    ATM_PRECIP_PASSES / ATM_PRECIP_RELAX. A new clean == old; C PASSES=30 RELAX=0.5 MUST differ.
# 2) 600 -> 620 restarts, -O2 build of the same source (cli/atm_conv), 6 x ${NT:-4} threads, diagnostics on every arm
#    (ATM_RAIN_PASS_DIAG, ATM_SNOW_DIAG, ATM_SR_DIAG, ATM_CWB_DIAG + CWB_BANDS):
#      default branch (output_sgzb_ctl):  conv_ru   UPWIND=1   (step 1: the converged budget)
#                                         conv_rd30 PASSES=30 RELAX=0.5, conv_rd60 PASSES=60 RELAX=0.3 (step 2)
#      closure branch (output_qh600):     conv_qu, conv_qd30, conv_qd60 (same three, qh600 setup)
# PRE-REGISTERED:
#  (a) the damped arms CONVERGE: |pass N - pass N-1| small against the band's rain, and d30 ~ d60 (if not, the shipped
#      formulation has no fixed point at all -- itself a result).
#  (b) the damped fixed point is compared with upwind: if both land at ~10 % of NASA, "the converged microphysics
#      rains a tenth of Earth" is robust; if the damped one is near the shipped ~988, upwind is the wrong form.
#  (c) step 1: from SNOW DIAG / SR DIAG on the converged arms, the rain-flux sources against S_ev per band -- the
#      recorded suspect is S_ev (demand 272 % of sources on 09-01); ATM_RAIN_AREA = 0.10 was FITTED on the
#      oscillating branch and is the first constant to re-examine.
set -u; cd "$(dirname "$0")"; rm -f CONV0928_DONE
. ./verify_lib.sh
[ -e ../cli/conv_old_atm ] || ./build_o0.sh 92f0dfa conv_old || { touch CONV0928_DONE; exit 1; }
[ -e ../cli/conv_new_atm ] || ./build_o0.sh HEAD    conv_new || { touch CONV0928_DONE; exit 1; }
for t in old new ctl; do mkdir output_vconv_$t || { touch CONV0928_DONE; exit 1; }; done
arm vconv_old ../cli/conv_old_atm config_vconv_old.xml
arm vconv_new ../cli/conv_new_atm config_vconv_new.xml
arm vconv_ctl ../cli/conv_new_atm config_vconv_ctl.xml ATM_PRECIP_PASSES=30 ATM_PRECIP_RELAX=0.5
wait_arms
cmp_dirs    vconv_new vconv_old "A  OFF BRANCH at -O0 (new clean vs 92f0dfa)"
want_differ vconv_ctl vconv_new "C  CONTROL PASSES=30 RELAX=0.5 -- MUST differ"
BAD=$(sed -n '1,/A  OFF BRANCH/p' run_conv_0928.out | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "C  CONTROL.*PASS" run_conv_0928.out; then
    echo "byte check did not pass -- restarts NOT started"; touch CONV0928_DONE; exit 1; fi
[ -e ../cli/atm_conv ] || { (cd .. && make atm > python/build_conv.log 2>&1) && cp ../cli/atm ../cli/atm_conv; }
for t in ru rd30 rd60 qu qd30 qd60; do mkdir output_conv_$t || { touch CONV0928_DONE; exit 1; }; done
for t in ru rd30 rd60; do cp output_sgzb_ctl/atm_restart_0Ma_600.bin output_conv_$t/; done
for t in qu qd30 qd60; do cp output_qh600/atm_restart_0Ma_600.bin    output_conv_$t/; done
QH="ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0"
D="ATM_RAIN_PASS_DIAG=1 ATM_SNOW_DIAG=1 ATM_SR_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1"
echo "restarts start $(date +%H:%M)"
run(){ ( env OMP_NUM_THREADS=${NT:-4} $D $2 ../cli/atm_conv config_conv_$1.xml > conv_$1.log 2>&1
         echo "conv_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' conv_$1.log)  $(date +%H:%M)" ) & }
run ru   "ATM_PRECIP_UPWIND=1"
run rd30 "ATM_PRECIP_PASSES=30 ATM_PRECIP_RELAX=0.5"
run rd60 "ATM_PRECIP_PASSES=60 ATM_PRECIP_RELAX=0.3"
run qu   "$QH ATM_PRECIP_UPWIND=1"
run qd30 "$QH ATM_PRECIP_PASSES=30 ATM_PRECIP_RELAX=0.5"
run qd60 "$QH ATM_PRECIP_PASSES=60 ATM_PRECIP_RELAX=0.3"
wait
for t in ru rd30 rd60 qu qd30 qd60; do echo "== conv_$t"; grep "RAIN PASS DIAG" conv_$t.log | tail -8
  grep "SNOW DIAG\] \(rain at ground\|rain-flux sources\|S_ev\)" conv_$t.log | tail -3; grep "SR DIAG" conv_$t.log | tail -3
  grep "by |latitude|" conv_$t.log | tail -1; grep "model .*NASA .*bias" conv_$t.log | tail -1 | cut -c1-150; done
touch CONV0928_DONE
