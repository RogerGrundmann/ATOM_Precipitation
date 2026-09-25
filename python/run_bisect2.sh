#!/bin/bash
# BISECTION of the def600 precipitation runaway (2026-09-25). 100 from scratch, cli/atm_def600 (= d76dd27,
# the binary def600 ran), 8 arms x 3 threads concurrent. Each arm reverts ONE 2026-09-24 default flip.
# def600 reference (24 threads): Precip 860 @40, 1401 @80, climbing; g_p = 0 at 600; CWB evaporation < 0.
#   bis_ctl  nothing reverted (same-thread-count control)
#   bis_wc   ATM_WATER_CLOSURE=0     bis_sf  ATM_EVAP_STRIDE_FIX=0    bis_sp  ATM_SEAM_PERIODIC=0
#   bis_ms   ATM_MC_S_NDIM=0         bis_hs  ATM_HYDRO_SPLIT=0        bis_rt  ATM_RAD_TOPO=0
#   bis_old  ALL of the above + ATM_RH_MIN_PTOP=490 -- the 09-23 default; must NOT run away, or the
#           cause is not among these flips.
# Score (score_bisect.py): Precip at 20/40/60/80/100 and its slope over 40-100, r, sigma, PW, bands,
# g_p / P_conv, CWB evaporation bucket sign.
set -u; cd "$(dirname "$0")"; rm -f BISECT2_DONE
echo "start $(date +%H:%M)"
run(){ local tag=$1; shift
       ( env OMP_NUM_THREADS=3 ATM_CWB_DIAG=1 ATM_MC_CAP_DIAG=1 "$@" ../cli/atm_def600 config_$tag.xml > $tag.log 2>&1
         echo "$tag exit $?  NaN $(grep -c 'NaN/Inf DETECTED' $tag.log)  $(date +%H:%M)" ) & }
run bis_ctl
run bis_wc  ATM_WATER_CLOSURE=0
run bis_sf  ATM_EVAP_STRIDE_FIX=0
run bis_sp  ATM_SEAM_PERIODIC=0
run bis_ms  ATM_MC_S_NDIM=0
run bis_hs  ATM_HYDRO_SPLIT=0
run bis_rt  ATM_RAD_TOPO=0
run bis_old ATM_WATER_CLOSURE=0 ATM_EVAP_STRIDE_FIX=0 ATM_SEAM_PERIODIC=0 ATM_MC_S_NDIM=0 ATM_HYDRO_SPLIT=0 ATM_RAD_TOPO=0 ATM_RH_MIN_PTOP=490
wait
touch BISECT2_DONE
