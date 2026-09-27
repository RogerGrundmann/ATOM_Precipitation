#!/bin/bash
# ATM_PRECIP_UPWIND -- the rain-iteration repair (2026-09-27). Both builds -O2.
# 1) byte check, 1 thread, nm = 20 from scratch: old = cli/atm_ret, new = cli/atm_upw2 (upwind knob + ZeroCat/OneCat/ThreeCat dropped).
#    A new clean == old;  C new + ATM_PRECIP_UPWIND=1 MUST differ.
# 2) waits for snw600 + ohs600 to free the machine, then 600 -> 640 restarts, 6 x ${NT:-4} threads, with
#    ATM_RAIN_PASS_DIAG + ATM_SNOW_DIAG + ATM_CWB_DIAG (+ CWB_BANDS) on every arm:
#    DEFAULT branch from output_sgzb_ctl:  upw_r0 (off), upw_r1 (UPWIND=1), upw_r13 (UPWIND=1 + SNOW_WINDOW=3)
#    CLOSURE branch from output_qh600:     upw_q0, upw_q1, upw_q13 (same three)
# PRE-REGISTERED:
#  (1) SELF-CHECK: in every UPWIND=1 arm, RAIN PASS DIAG |pass 3 - pass 2| = 0.000 exactly in all four bands.
#  (2) the question: where does the converged rain land? Shipped pass 3 gives 0-15 ~3465 / 35-65 ~201 mm/a,
#      pass 1 ~370 / ~141. If the converged 0-15 value is far below 3465, the tropical spike (2.3x NASA) was the
#      oscillation, not the circulation. No expectation registered for the value itself.
#  (3) r13 vs r1: with the iteration converged, does keeping and melting the snow ADD rain now (the claim the
#      snow window was written on), or does it still take rain away?
#  (4) exit 0, zero NaN; floor injection (printed every run) stays ~0; default-branch arms (r*) re-pin the surface
#      humidity, so the closure arms (q*) carry the water-conserving amounts.
set -u; cd "$(dirname "$0")"; rm -f UPW40_DONE
. ./verify_lib.sh
for t in old new ctl; do mkdir output_vupw_$t || { touch UPW40_DONE; exit 1; }; done
arm vupw_old ../cli/atm_ret config_vupw_old.xml
arm vupw_new ../cli/atm_upw2 config_vupw_new.xml
arm vupw_ctl ../cli/atm_upw2 config_vupw_ctl.xml ATM_PRECIP_UPWIND=1
wait_arms
cmp_dirs    vupw_new vupw_old "A  OFF BRANCH (new clean vs old, both -O2)"
want_differ vupw_ctl vupw_new "C  CONTROL PRECIP_UPWIND=1 -- MUST differ"
BAD=$(sed -n '1,/A  OFF BRANCH/p' run_upw40.out | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "C  CONTROL.*PASS" run_upw40.out; then
    echo "byte check did not pass -- restarts NOT started"; touch UPW40_DONE; exit 1; fi
until [ -f SNW600_DONE ] && [ -f OHS600_DONE ]; do sleep 30; done
for t in r0 r1 r13 q0 q1 q13; do mkdir output_upw_$t || { touch UPW40_DONE; exit 1; }; done
for t in r0 r1 r13; do cp output_sgzb_ctl/atm_restart_0Ma_600.bin output_upw_$t/; done
for t in q0 q1 q13; do cp output_qh600/atm_restart_0Ma_600.bin    output_upw_$t/; done
QH="ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0"
echo "restarts start $(date +%H:%M)"
run(){ ( env OMP_NUM_THREADS=${NT:-4} ATM_RAIN_PASS_DIAG=1 ATM_SNOW_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 $2 \
             ../cli/atm_upw2 config_upw_$1.xml > upw_$1.log 2>&1
         echo "upw_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' upw_$1.log)  $(date +%H:%M)" ) & }
run r0 ""; run r1 "ATM_PRECIP_UPWIND=1"; run r13 "ATM_PRECIP_UPWIND=1 ATM_SNOW_WINDOW=3"
run q0 "$QH"; run q1 "$QH ATM_PRECIP_UPWIND=1"; run q13 "$QH ATM_PRECIP_UPWIND=1 ATM_SNOW_WINDOW=3"
wait
for t in r0 r1 r13 q0 q1 q13; do echo "== upw_$t"; grep "RAIN PASS DIAG" upw_$t.log | tail -6; grep "by |latitude|" upw_$t.log | tail -1
  grep "model .*NASA .*bias" upw_$t.log | tail -1 | cut -c1-150; grep -i "floor inject" upw_$t.log | tail -1 | cut -c1-90; done
touch UPW40_DONE
