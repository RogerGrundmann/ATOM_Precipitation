#!/bin/bash
# 2026-09-28: ATM_SEAM_Q_DIAG (print-only) -- where does the seam water leak ATM_SEAM_Q_CONSERVE leaves (fx4_sqcon:
# 1317 -> 959 mm/a, 3052 in 15-35 deg) come from? Starts when run_rp0928.sh is done.
# 1) -O0 byte check (both sides -O0): old = cli/conv_new_atm (HEAD c1aae23, built this morning), new = WORKTREE.
#    A  new clean == old;  B  new + DIAG == new (print-only);  C  new CONSERVE=1 + DIAG == new CONSERVE=1;
#    D  the diag must actually print (grep), else the checks above are vacuous.
# 2) 600 -> 620 restarts, -O2 build of the worktree (cli/atm_sqd), 3 x 4 threads, CWB diag on:
#      sqd_qoff  qh600 setup (closure, moisture filter off), DIAG
#      sqd_qon   same + ATM_SEAM_Q_CONSERVE=1
#      sqd_def   default branch from sgzb_ctl, DIAG
# PRE-REGISTERED: in sqd_qon the conserved set (A) sums to ~0 and the leak is in B (seam air, land neighbour),
# concentrated in 15-35 deg (0 deg longitude crosses the Sahara). If A carries it, the conserving step itself is
# wrong (mass mismatch: it uses the post-average r_humid). The A* clip row sizes the no-negative clip.
set -u; cd "$(dirname "$0")"; rm -f SQD0928_DONE
until [ -e RP0928_DONE ]; do sleep 60; done
. ./verify_lib.sh
[ -e ../cli/sqd_new_atm ] || ./build_o0.sh WORKTREE sqd_new || { touch SQD0928_DONE; exit 1; }
for t in old new d c cd; do mkdir output_vsqd_$t || { touch SQD0928_DONE; exit 1; }
  sed "s#output_vconv_new/#output_vsqd_$t/#" config_vconv_new.xml > config_vsqd_$t.xml; done
arm vsqd_old ../cli/conv_new_atm config_vsqd_old.xml
arm vsqd_new ../cli/sqd_new_atm  config_vsqd_new.xml
arm vsqd_d   ../cli/sqd_new_atm  config_vsqd_d.xml   ATM_SEAM_Q_DIAG=1
arm vsqd_c   ../cli/sqd_new_atm  config_vsqd_c.xml   ATM_SEAM_Q_CONSERVE=1
arm vsqd_cd  ../cli/sqd_new_atm  config_vsqd_cd.xml  ATM_SEAM_Q_CONSERVE=1 ATM_SEAM_Q_DIAG=1
wait_arms
for p in "old:new" "d:new" "cd:c"; do a=${p%%:*}; b=${p##*:}
  for f in output_vsqd_$a/RUN_CONFIG.txt output_vsqd_$b/RUN_CONFIG.txt; do sed -i 's#output_vsqd_[a-z]*/#OUT/#' $f; done; done
cmp_dirs vsqd_new vsqd_old "A  OFF BRANCH at -O0 (worktree clean vs HEAD)"
cmp_dirs vsqd_d   vsqd_new "B  DIAG is print-only"
cmp_dirs vsqd_cd  vsqd_c   "C  DIAG is print-only under CONSERVE=1"
N=$(grep -c 'SEAM Q DIAG\] total' vsqd_d.log); echo "D  diag lines printed: $N  -> $([ "$N" -gt 0 ] && echo PASS || echo FAIL)"
if grep -E '^(A|B|C) .*FAIL' run_sqd0928.out >/dev/null || [ "$N" = 0 ]; then
    echo "byte check did not pass -- restarts NOT started"; touch SQD0928_DONE; exit 1; fi
(cd .. && make atm > python/build_sqd2.log 2>&1) && cp -n ../cli/atm ../cli/atm_sqd || { touch SQD0928_DONE; exit 1; }
for t in qoff qon def; do mkdir output_sqd_$t || { touch SQD0928_DONE; exit 1; }
  sed "s#output_conv_qu/#output_sqd_$t/#" config_conv_qu.xml > config_sqd_$t.xml; done
cp output_qh600/atm_restart_0Ma_600.bin output_sqd_qoff/; cp output_qh600/atm_restart_0Ma_600.bin output_sqd_qon/
cp output_sgzb_ctl/atm_restart_0Ma_600.bin output_sqd_def/
QH="ATM_WATER_CLOSURE=1 ATM_DAMP_Q_VERT=0 ATM_DAMP_Q_HORIZ=0 ATM_DAMP_T_VERT=0"
D="ATM_SEAM_Q_DIAG=1 ATM_CWB_DIAG=1 ATM_CWB_BANDS=1"
echo "restarts start $(date +%H:%M)  atm_sqd $(md5sum < ../cli/atm_sqd | cut -c1-8)"
run(){ ( env OMP_NUM_THREADS=${NT:-4} $D $2 ../cli/atm_sqd config_sqd_$1.xml > sqd_$1.log 2>&1
         echo "sqd_$1 exit $?  NaN $(grep -c 'NaN/Inf DETECTED' sqd_$1.log)  $(date +%H:%M)" ) & }
run qoff "$QH"
run qon  "$QH ATM_SEAM_Q_CONSERVE=1"
run def  ""
wait
for t in qoff qon def; do echo "== sqd_$t"; grep 'SEAM Q DIAG' sqd_$t.log | tail -7; grep 'CWB\] BC:phi(seam)' sqd_$t.log | tail -1; grep 'CWB-BANDS\] BC:phi(seam)' sqd_$t.log | tail -1; done
touch SQD0928_DONE
