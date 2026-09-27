#!/bin/bash
# THE SNOW WINDOW, SIZED AND REPAIRED (2026-09-27, prepared for the storm-track item; NOT QUEUED -- launch by hand).
# 600 -> 640 restarts, cli/atm_snw, 6 arms x ${NT:-4} threads concurrent, ATM_SNOW_DIAG=1 on all, gated on the byte check.
#   snw_r0/r1/r2/r3  DEFAULT branch from output_sgzb_ctl's checkpoint, SNOW_WINDOW = 0 / 1 / 2 / 3
#   snw_q0/q3        CLOSURE branch (qh600 setup) from output_qh600's checkpoint, SNOW_WINDOW = 0 / 3
# WHY: TwoCat zeroes the snow flux outside -20..0 C. Snow reaching the melting level is DELETED (S_s_melt, S_shed
# read the same-level flux the window set to 0: structurally zero), and snow below -20 C likewise, while S_i_au has
# no temperature gate. Storm-track precipitation on Earth is mostly ice-phase aloft that melts on the way down.
# PRE-REGISTERED:
#  (1) r0/q0 size it: "snow DELETED warm" in 35-65 deg. If it is < ~50 mm/a (against a 890 mm/a gap to NASA's 981)
#      the window is REFUTED as the storm-track lever; if it is hundreds of mm/a it is the lever.
#  (2) r1: "deleted warm" -> exactly 0 (self-check), S_s_melt > 0, and 35-65 rain at the ground rises by roughly
#      the deleted amount less what S_ev takes. r2: "deleted cold" -> exactly 0. r3: both.
#  (3) default-branch arms (r*) re-pin the surface humidity, so read them for MECHANISM; the closure arms (q*)
#      conserve water, so q3 - q0 is the amount that counts.
#  (4) exit 0, zero NaN; P_max cap: max precipitation cells.
set -u; cd "$(dirname "$0")"; rm -f SNW40_DONE
until [ -f VSNW_VERIFY_DONE ]; do sleep 30; done
BAD=$(sed -n '1,/A  OFF BRANCH/p' run_verify_vsnw.out | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BAD" != 0 ] || ! grep -q "A  OFF BRANCH" run_verify_vsnw.out || ! grep -q "D  SNOW_DIAG.*PASS" run_verify_vsnw.out \
   || [ "$(grep -c 'CONTROL.*PASS' run_verify_vsnw.out)" != 2 ]; then
    echo "byte check did not pass -- NOT started"; cat run_verify_vsnw.out; touch SNW40_DONE; exit 1; fi
for t in r0 r1 r2 r3 q0 q3; do mkdir output_snw_$t || { touch SNW40_DONE; exit 1; }; done
for t in r0 r1 r2 r3; do cp output_sgzb_ctl/atm_restart_0Ma_600.bin output_snw_$t/; done
for t in q0 q3;       do cp output_qh600/atm_restart_0Ma_600.bin    output_snw_$t/; done
QH="ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0"
echo "start $(date +%H:%M)"
run(){ ( env OMP_NUM_THREADS=${NT:-4} ATM_SNOW_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 $2 ../cli/atm_snw config_snw_$1.xml > snw_$1.log 2>&1
         echo "snw_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' snw_$1.log)  $(date +%H:%M)" ) & }
run r0 "ATM_SNOW_WINDOW=0"; run r1 "ATM_SNOW_WINDOW=1"; run r2 "ATM_SNOW_WINDOW=2"; run r3 "ATM_SNOW_WINDOW=3"
run q0 "$QH ATM_SNOW_WINDOW=0"; run q3 "$QH ATM_SNOW_WINDOW=3"
wait
for t in r0 r1 r2 r3 q0 q3; do echo "== snw_$t"; grep "SNOW DIAG" snw_$t.log | tail -9; grep "by |latitude|" snw_$t.log | tail -1; done
touch SNW40_DONE
