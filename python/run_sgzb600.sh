#!/bin/bash
# B.10b arms on the CORRECTED default (ATM_WATER_CLOSURE back OFF, 2026-09-25): ATM_MC_SGZ + ATM_MC_ENTR.
# 600 from scratch, binary cli/atm_wcr (md5 314e7e99), 4 arms x 6 threads concurrent. ATM_CWB_DIAG +
# ATM_MC_CAP_DIAG on. Gated on run_verify_wcr.sh (closure default revert byte check).
#   sgzb_ctl   the corrected default (no 600 run of it exists yet -- this is also that run)
#   sgzb_e     ATM_MC_ENTR=1e-4                      entrainment alone
#   sgzb_sg    ATM_MC_SGZ=1                          g*z alone -- the predicted sign flip
#   sgzb_sge   ATM_MC_SGZ=1 ATM_MC_ENTR=1e-4         the pair
# PRE-REGISTERED (parcel on the 09-24 default: 1.6-3.0 K warm, mixing fraction D_u*dz/M_u 0.66-1.67; on
# the closure-off branch bis_wc generated g_p ~1058 mm/a at iter 100, P_conv 0.22):
#   0. ctl: no runaway -- precip flattens (bis_wc: 982 @100, slope 3.6 and falling); bands near 09-22.
#   1. e: mixing fraction ~20x smaller; with no g*z the parcel keeps ~cloud-base cp*T and runs FAR too
#      warm aloft (tens of K) -- s_u - s grows with height.
#   2. sg: s_u - s NEGATIVE through the cloud layer -> MC_t flux half flips sign (cools aloft).
#   3. sge: parcel follows a (moist) adiabat, a few K warm; g_p changes; does P_conv leave ~0.2 mm/a?
#   4. every arm: exit 0, zero NaN through 155/357/483; MC_t truncation; max u/v/w vs ctl; precip
#      global/bands/r/sigma vs ctl at 600 (convection is tropical: 0-15 is where a change should show).
# Instrument: parcel vs environment from the zonal-87 slice at 100/.../600 (s_u, s, M_u, D_u).
set -u; cd "$(dirname "$0")"; rm -f SGZB600_DONE
until [ -f WCR_VERIFY_DONE ]; do sleep 30; done
# every DIFFERS line printed before the C verdict must be RUN_CONFIG.txt, and A/B must each report 13/1 or PASS
BADALL=$(sed -n '1,/C  CONTROL/p' run_verify_wcr.out | sed '/C  CONTROL/,$d' | grep 'DIFFERS:' | grep -vc 'RUN_CONFIG.txt')
if [ "$BADALL" != 0 ] || ! grep -q "C  CONTROL.*PASS" run_verify_wcr.out \
   || [ "$(grep -cE '^(A|B)  .*(-> PASS|identical 13 / differing 1)' run_verify_wcr.out)" != 2 ]; then
    echo "byte check did not pass -- arms NOT started"; cat run_verify_wcr.out; touch SGZB600_DONE; exit 1; fi
for t in ctl e sg sge; do mkdir output_sgzb_$t || { echo "output_sgzb_$t exists -- NOT started"; touch SGZB600_DONE; exit 1; }; done
echo "start $(date +%H:%M)"
run(){ local tag=$1; shift
       ( env OMP_NUM_THREADS=6 ATM_CWB_DIAG=1 ATM_MC_CAP_DIAG=1 "$@" ../cli/atm_wcr config_$tag.xml > $tag.log 2>&1
         echo "$tag exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $tag.log)  $(date +%H:%M)" ) & }
run sgzb_ctl
run sgzb_e   ATM_MC_ENTR=1.0e-4
run sgzb_sg  ATM_MC_SGZ=1
run sgzb_sge ATM_MC_SGZ=1 ATM_MC_ENTR=1.0e-4
wait
touch SGZB600_DONE
