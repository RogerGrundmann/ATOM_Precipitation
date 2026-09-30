#!/bin/bash
# 2026-09-30: bisect the batch-1 plan C signed-zero S_r difference. Each variant = HEAD + a subset of the batch-1
# diff (CloudFraction.h and Knobs.h never included, so every subset compiles), -O0 atm only, nm 1, 1 thread,
# compared at iteration 0 against output_vretb1_b (the old binary).
set -u; cd "$(dirname "$0")/.."
declare -A V=( [v1]="atmosphere/TwoCatIceScheme.h" [v2]="atmosphere/SaturationAdjustment.h"
               [v3]="atmosphere/InitValues_Atm.cpp"
               [v4]="atmosphere/MoistConvection.h atmosphere/MultiLayerRadiation.h atmosphere/UtilsAtm.h atmosphere/ThermoAtm.h" )
build(){ local n=$1; shift; local W; W=$(mktemp -d "${TMPDIR:-/tmp}/atom_bis_XXXXXX")
  git worktree add --detach "$W" HEAD >/dev/null && git diff HEAD -- "$@" > "$W/.p" && (cd "$W" && git apply .p && make -j4 OPT=-O0 atm > b.log 2>&1) \
    && cp "$W/cli/atm" cli/bis_$n; git worktree remove --force "$W"; }
for n in v1 v2 v3 v4; do build $n ${V[$n]} & done; wait; git worktree prune
cd python
for n in v1 v2 v3 v4; do [ -e ../cli/bis_$n ] || { echo "$n build failed"; continue; }
  mkdir output_bis_$n || { echo "output_bis_$n exists -- abort"; exit 1; }; sed -e "s#output_vconv_new/#output_bis_$n/#" -e 's#<nm>20</nm>#<nm>1</nm>#' config_vconv_new.xml > config_bis_$n.xml
  ( env OMP_NUM_THREADS=1 ../cli/bis_$n config_bis_$n.xml > bis_$n.log 2>&1
    for f in 0Ma_smooth_Atm_radial_0_0.vtk 0Ma_smooth_Atm_zonal_87_0.vtk 0Ma_smooth_Transfer_Atm_0.vwtp; do
      cmp -s output_bis_$n/$f output_vretb1_b/$f && r=same || r=DIFF; echo "$n ${V[$n]%% *}...: $f $r"; done ) &
done; wait; touch BISECT_B1_DONE
