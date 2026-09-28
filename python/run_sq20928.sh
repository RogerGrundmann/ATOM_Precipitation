#!/bin/bash
# 2026-09-28: ATM_SEAM_Q_CONSERVE=2 -- also conserve in the seam cells with a LAND neighbour (the 72 cells that carry
# 99.6 % of the leak =1 leaves, ATM_SEAM_Q_DIAG / run_sqd0928.sh).
# 1) -O0 byte check, 1 thread, nm 20 from scratch: old = HEAD (caf17b7) -O0, new = WORKTREE -O0
#      A  new clean == old clean;  B  new CONSERVE=1 == old CONSERVE=1 (mode 1 untouched)
#      C  CONTROL new CONSERVE=2 != new CONSERVE=1
# 2) -O2 (cli/atm_sq2, built from the worktree), 600 -> 620 restarts, 4 threads, SEAM_Q_DIAG + CWB:
#      sq2_qon  qh600 setup + CONSERVE=2    against sqd_qoff / sqd_qon (atm_sqd; differs only by the null TURB_SIN flip)
#      sq2_def  default branch + CONSERVE=2 against sqd_def
# PRE-REGISTERED: in sq2_qon category B drops ~100x like A did, the CWB BC:phi(seam) row falls from ~980 to
# ~0-20 mm/a and 15-35 deg from ~3100 to ~0; the no-negative clip is the residual. sq2_def: the -2.27e4 mm/a
# default-branch leak falls to at most the A-category part.
set -u; cd "$(dirname "$0")"; rm -f SQ20928_DONE
. ./verify_lib.sh
[ -e ../cli/sq2_old_atm ] || ./build_o0.sh HEAD     sq2_old || { touch SQ20928_DONE; exit 1; }
[ -e ../cli/sq2_new_atm ] || ./build_o0.sh WORKTREE sq2_new || { touch SQ20928_DONE; exit 1; }
for t in o0 n0 o1 n1 n2; do mkdir output_vsq2_$t || { touch SQ20928_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vsq2_$t/#" config_vconv_new.xml > config_vsq2_$t.xml; done
arm vsq2_o0 ../cli/sq2_old_atm config_vsq2_o0.xml
arm vsq2_n0 ../cli/sq2_new_atm config_vsq2_n0.xml
arm vsq2_o1 ../cli/sq2_old_atm config_vsq2_o1.xml ATM_SEAM_Q_CONSERVE=1
arm vsq2_n1 ../cli/sq2_new_atm config_vsq2_n1.xml ATM_SEAM_Q_CONSERVE=1
arm vsq2_n2 ../cli/sq2_new_atm config_vsq2_n2.xml ATM_SEAM_Q_CONSERVE=2
wait_arms
for t in o0 n0 o1 n1 n2; do sed -i 's#output_vsq2_[a-z0-9]*/#OUT/#' output_vsq2_$t/RUN_CONFIG.txt; done
cmp_dirs vsq2_n0 vsq2_o0 "A  OFF BRANCH at -O0"
cmp_dirs vsq2_n1 vsq2_o1 "B  MODE 1 UNCHANGED"
want_differ vsq2_n2 vsq2_n1 "C  CONTROL CONSERVE=2 vs =1 -- MUST differ"
if grep -E '^(A|B) .*FAIL' run_sq20928.out >/dev/null || ! grep -q 'C  CONTROL.*PASS' run_sq20928.out; then
    echo "byte check did not pass -- restarts NOT started"; touch SQ20928_DONE; exit 1; fi
for t in qon def; do mkdir output_sq2_$t || { touch SQ20928_DONE; exit 1; }
  sed "s#output_conv_qu/#output_sq2_$t/#" config_conv_qu.xml > config_sq2_$t.xml; done
cp output_qh600/atm_restart_0Ma_600.bin output_sq2_qon/; cp output_sgzb_ctl/atm_restart_0Ma_600.bin output_sq2_def/
QH="ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0"
D="ATM_SEAM_Q_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1 ATM_SEAM_Q_CONSERVE=2"
echo "restarts start $(date +%H:%M)  atm_sq2 $(md5sum < ../cli/atm_sq2 | cut -c1-8)"
run(){ ( env OMP_NUM_THREADS=${NT:-4} $D $2 ../cli/atm_sq2 config_sq2_$1.xml > sq2_$1.log 2>&1
         echo "sq2_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' sq2_$1.log)  $(date +%H:%M)" ) & }
run qon "$QH"
run def ""
wait
for t in qon def; do echo "== sq2_$t"; grep 'SEAM Q DIAG' sq2_$t.log | tail -7; grep 'CWB\] BC:phi(seam)' sq2_$t.log | tail -1; grep 'CWB-BANDS\] BC:phi(seam)' sq2_$t.log | tail -1
  grep "by |latitude|" sq2_$t.log | tail -1; done
touch SQ20928_DONE
