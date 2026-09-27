#!/bin/bash
# ATM_PRECIP_UPWIND restarts (2026-09-27), the second half of run_upw40.sh without its -O2 byte check, which the
# -O0 check (run_verify_o0.sh: 13/13 identical) superseded. Binary cli/atm_92f = HEAD 92f0dfa at -O2. Same arms and
# pre-registration as run_upw40.sh.
set -u; cd "$(dirname "$0")"; rm -f UPW40_DONE
for t in r0 r1 r13 q0 q1 q13; do mkdir output_upw_$t || { touch UPW40_DONE; exit 1; }; done
for t in r0 r1 r13; do cp output_sgzb_ctl/atm_restart_0Ma_600.bin output_upw_$t/; done
for t in q0 q1 q13; do cp output_qh600/atm_restart_0Ma_600.bin    output_upw_$t/; done
QH="ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0"
echo "restarts start $(date +%H:%M)"
run(){ ( env OMP_NUM_THREADS=${NT:-4} ATM_RAIN_PASS_DIAG=1 ATM_SNOW_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 $2 \
             ../cli/atm_92f config_upw_$1.xml > upw_$1.log 2>&1
         echo "upw_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' upw_$1.log)  $(date +%H:%M)" ) & }
run r0 ""; run r1 "ATM_PRECIP_UPWIND=1"; run r13 "ATM_PRECIP_UPWIND=1 ATM_SNOW_WINDOW=3"
run q0 "$QH"; run q1 "$QH ATM_PRECIP_UPWIND=1"; run q13 "$QH ATM_PRECIP_UPWIND=1 ATM_SNOW_WINDOW=3"
wait
for t in r0 r1 r13 q0 q1 q13; do echo "== upw_$t"; grep "RAIN PASS DIAG" upw_$t.log | tail -6; grep "by |latitude|" upw_$t.log | tail -1
  grep "model .*NASA .*bias" upw_$t.log | tail -1 | cut -c1-150; grep -i "floor inject" upw_$t.log | tail -1 | cut -c1-90; done
touch UPW40_DONE
