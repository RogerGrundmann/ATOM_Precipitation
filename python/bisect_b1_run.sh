#!/bin/bash
# run the four bisection binaries (cli/bis_v1..v4, built by bisect_b1.sh) -- nm 1, 1 thread, compare iteration 0
set -u; cd "$(dirname "$0")"; rm -f BISECT_B1_DONE
declare -A V=( [v1]="TwoCat" [v2]="SaturationAdjustment" [v3]="InitValues_Atm" [v4]="MoistConv+MLR+UtilsAtm+ThermoAtm" )
for n in v1 v2 v3 v4; do mkdir output_bisb1_$n || { echo "output_bisb1_$n exists -- abort"; touch BISECT_B1_DONE; exit 1; }
  sed -e "s#output_vconv_new/#output_bisb1_$n/#" -e 's#<nm>20</nm>#<nm>1</nm>#' config_vconv_new.xml > config_bisb1_$n.xml; done
for n in v1 v2 v3 v4; do
  ( env OMP_NUM_THREADS=1 ../cli/bis_$n config_bisb1_$n.xml > bisb1_$n.log 2>&1
    for f in 0Ma_smooth_Atm_radial_0_0.vtk 0Ma_smooth_Atm_zonal_87_0.vtk 0Ma_smooth_Transfer_Atm_0.vwtp; do
      cmp -s output_bisb1_$n/$f output_vretb1_b/$f && r=same || r=DIFF; echo "$n ${V[$n]}: $f $r"; done ) &
done; wait; touch BISECT_B1_DONE
