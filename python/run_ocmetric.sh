#!/bin/bash
# THE OCEAN HORIZONTAL-METRIC REPAIR, FROM SCRATCH, AS A FOUR-STEP LADDER.
#
# THE DEFECT. `RungeKutta_Hyd_Turb` sets geo.rm = rad.z[i], rad.z runs 1.000-2.000 and one unit
# is L_hyd = 200 m, so the implied PLANETARY radius in every horizontal derivative is 200-400 m.
# Measured on the model's own field: the code's |horizontal advection|/|Coriolis| is 18 972x the
# same ratio formed physically. The atmosphere does not have this defect -- ATM_METRIC_RADIUS
# DEFAULTS TO r_Earth -- and the mechanism is already ported to the ocean verbatim
# (metricRadius(rm) = m_metric_r0 + (rm - rad.z[0]), identity when off). **So nothing has to be
# WRITTEN. The asymmetry between the two models is one default.**
#
# ⚠⚠ FROM SCRATCH IS MANDATORY AND HAS NEVER BEEN DONE. Every recorded metric arm (r0/r1/r2,
# m0, w0/w1/w2, kr0/kr1, va0/va1/va2, f0/f1/s1/s2, von/voff) is a RESTART from
# `hyd_restart_0Ma_300.bin`, and that checkpoint PREDATES `HYD_SSS_FILL` (default ON since
# 2026-09-08). It therefore carries the SSS sentinel field: 13.73 % of fluid cells below 5 psu,
# ocean mean salinity 28.6 psu against Earth's 34.7. `load_state` overwrites `c`, so no restart
# can see the repair. Consequences for this ladder, both of which are the point of running it:
#   * the recorded "2.3x grid-scale noise" and "bottom-intensified profile" were measured on an
#     ocean that was 6 psu fresh in the mean, and have never been checked on a correct one;
#   * `HYD_PHYDRO_SALT` could not have worked there anyway -- its plausibility floor was
#     rejecting 12.35 % of cells and falling back to `r_water`, which MANUFACTURES the very
#     horizontal density gradient the field is being repaired to provide. With SSS_FILL on the
#     floor rejects 93 % fewer cells, so step 4 below is being given its prerequisite for the
#     first time.
#
# THE LADDER -- one change per rung, so every effect is attributable:
#   oc_ctl   shipped everything                    the reference PROFILE and the noise floor
#   oc_met   + HYD_METRIC_RADIUS=6370 HYD_RUN_NEUMANN=1
#   oc_vis   + HYD_A_H_BIHARM=3.0e18
#   oc_full  + HYD_PHYDRO_SALT=1 HYD_BAROCLINIC_PGF=1.0
#
# WHY RUN_NEUMANN IS NOT OPTIONAL ON RUNG 2. The metric alone gives
# "NaN/Inf DETECTED stage=post-RK4 iter=1 field=u". That was first diagnosed as solver
# conditioning and THAT DIAGNOSIS IS RETRACTED: `PressureSolverHyd::run()`'s radial pressure BC
# is a CUBIC EXTRAPOLATION, p[0] = p[3] - 3p[2] + 3p[1] -- a zero-third-derivative SMOOTHNESS
# condition where a projection requires dp/dn = 0 -- and it permits flow through the seafloor
# and the lid. It is an exact NULL at the shipped metric (max|u| 0.016243 both ways), which is
# both the attribution and the reason it went unfound for so long: at a 200 m metric radius the
# horizontal terms swamp the radial BC error by 2e4.
#
# WHY B = 3.0e18 m^4/s, DERIVED RATHER THAN CHOSEN. The metric error was SUPPLYING the ocean's
# horizontal viscosity -- implied 1.16e+06 m^2/s at the surface and 4.64e+06 at the seafloor,
# against 1e2..1e4 in a 1-degree ocean model. Correcting the metric removes it (inv_rm^2 ~ 4e8)
# and leaves nothing to damp grid-scale structure. The measured 2026-09-05 arm used
# B = 1.882e+17 (2dx e-folding 100 iterations) and cut the noise index only 9.0 %, because that
# index is NOT a pure 2dx mode -- the Laplacian at equal nominal strength got 28.7 % precisely
# because it also damps 3dx and 4dx, where k^4 has fallen off by 5x and 16x. So the biharmonic
# needs ~16x to reach 4dx, AND IT CAN AFFORD IT: at 1.882e+17 the domain-scale e-folding is
# 6313 days, so 16x still leaves ~396 days = 4.1e+08 iterations, i.e. no circulation cost over a
# run of 1000. PREDICTION TO CHECK AGAINST THE MODEL'S OWN [SCALES] PRINT: 2dx e-folding ~6.3
# iterations, domain scale ~396 d. Explicit stability wants B <= 2.3e+20, so this is 77x under.
# ⚠ Quote it as a NUMERICAL knob. A physical 1e4 m^2/s needs 1.5e6 iterations to act.
#
# ⚠ HYD_BAROCLINIC_PGF=1.0 IS ~5x THE CORIOLIS FORCE IT EXISTS TO BALANCE on the salt-aware
# field (median ratio 7.96 all cells, 5.13 on salty columns). Run it at 1.0 anyway: the strength
# sweep is EXACTLY LINEAR with identical stability at 0.1/0.2/0.5/1.0, and the excess is the
# FLOW's -- the density field implies a geostrophic current 2.4x what the model has -- so
# scaling the force down would be fitting a derived term to a wrong flow.
#
# FORCING: output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp, copied into each arm's directory,
# atm_transfer_iter = 600. Deliberately the SAME file every recorded ocean arm used, so the one
# axis still comparable to the record -- the shipped-metric control -- stays comparable.
# ⚠ It is a PRE-2026-09-12 atmosphere. And the file must be present: the HYD_T_FREEZE_SFC null
# of 2026-09-05 was an output directory seeded with the .bin alone, falling back to the NASA
# surface temperature, with the WARNING forty lines above the line being read.
#
# PRE-REGISTERED EXPECTATIONS:
#  1. oc_met reproduces the recorded failure FROM SCRATCH: rms u collapses ~2000x (that part is
#     the repair working), noise index up ~2x, and the velocity profile goes flat-to-bottom-
#     intensified. IF IT DOES NOT, the salinity field was carrying that result and every
#     recorded metric arm was measured on a defective ocean.
#  2. oc_vis controls the noise and costs no circulation. If the noise is still uncontrolled at
#     16x the previously measured strength, the index is not a grid-scale mode at all.
#  3. oc_full is the only rung that can fix the PROFILE. The bottom-intensification mechanism --
#     that the spurious radial velocity was doing the vertical momentum transport, and removing
#     it leaves nothing to shear the column -- is recorded as NOT MEASURED. This tests it. If the
#     profile is still inverted with a working baroclinic PGF, that mechanism is wrong.
#  4. STABILITY: the documented from-scratch hazard is the barotropic cold-start CFL runaway at
#     iterations 65-69, which is what the 120-iteration ramp exists for. Watch there, not at 155/
#     357/483 -- those are the ATMOSPHERE's failure points.
#  5. SCORE ON THE PROFILE, NOT ON max. `max` is one cell: the 2026-09-04 "ACC jet appears at
#     1.12 m/s" was a single noisy cell read as a feature and had to be retracted, while p50/p90/
#     p99 were identical between the arms. Fixed population of full-depth columns per level.
#
# 4 arms x 6 threads, concurrent. Waits for the PTOP sweep to release the cores.
set -u; cd "$(dirname "$0")"
rm -f OCMETRIC_DONE
echo "=== waiting for PTOP600_DONE $(date +%H:%M:%S)"
until [ -f PTOP600_DONE ]; do sleep 60; done
echo "=== cores free, launching ocean ladder $(date +%H:%M:%S)"
export OMP_NUM_THREADS=6
XFER=output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp
for t in oc_ctl oc_met oc_vis oc_full; do
  rm -rf output_$t; mkdir -p output_$t
  cp "$XFER" "output_$t/" || { echo "MISSING $XFER"; exit 1; }
done
( env                                                                                   ../cli/hyd_oc config_oc_ctl.xml  > oc_ctl.log  2>&1; echo "  oc_ctl  exit $? nan $(grep -ci 'nan' oc_ctl.log)  $(date +%H:%M:%S)" ) &
( env HYD_METRIC_RADIUS=6370 HYD_RUN_NEUMANN=1                                          ../cli/hyd_oc config_oc_met.xml  > oc_met.log  2>&1; echo "  oc_met  exit $? nan $(grep -ci 'nan' oc_met.log)  $(date +%H:%M:%S)" ) &
( env HYD_METRIC_RADIUS=6370 HYD_RUN_NEUMANN=1 HYD_A_H_BIHARM=3.0e18                    ../cli/hyd_oc config_oc_vis.xml  > oc_vis.log  2>&1; echo "  oc_vis  exit $? nan $(grep -ci 'nan' oc_vis.log)  $(date +%H:%M:%S)" ) &
( env HYD_METRIC_RADIUS=6370 HYD_RUN_NEUMANN=1 HYD_A_H_BIHARM=3.0e18 HYD_PHYDRO_SALT=1 HYD_BAROCLINIC_PGF=1.0 ../cli/hyd_oc config_oc_full.xml > oc_full.log 2>&1; echo "  oc_full exit $? nan $(grep -ci 'nan' oc_full.log)  $(date +%H:%M:%S)" ) &
wait
echo "=== ocean ladder done $(date +%H:%M:%S)"
touch OCMETRIC_DONE
