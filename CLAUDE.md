# ATOM_Precipitation — modern Earth, and the tree the family forked from

Coupled atmosphere + ocean on a spherical shell: finite-difference Navier-Stokes, RK4,
exponentially stretched vertical coordinate, optional k-eps / k-omega / k-omega-SST closures,
COSMO ice schemes and a saturation adjustment. Paleo time slices driven by Scotese
reconstructions. **This is the PARENT of ATHAD, ATHAD_COND and ATHAD_PERID** — class, file and
function names are deliberately identical across those forks so fixes cherry-pick in both
directions.

## Build, run

```bash
make atm            # -> cli/atm          make hyd -> cli/hyd
cd python && OMP_NUM_THREADS=8 ../cli/atm config_atm.xml
```

`make` regenerates the parameter bindings from `param.py` first. **Every object depends on the
`Makefile`**, deliberately — see the build hazard below.

## READ THIS BEFORE QUOTING ANY NUMBER FROM THIS TREE

**EVERY FIGURE RECORDED BEFORE 2026-08-26 CAME FROM A BINARY WITH THREE DATA RACES IN IT**, and
is not reproducible. Nothing is retracted — the races moved last digits — but the CURE moves
more than the races did, because red-black relaxation is a different sweep order from
lexicographic Gauss-Seidel: `residuum_atm` shifts 2.3e-6 while **`max u-component` moves
0.89 %** at four iterations. A figure quoted to better than ~1 % from before that date should be
re-measured rather than cited. See README, *Correctness fixes ported from ATHAD*.

**AND EVERY `Psi` NUMBER PREDATING THAT DATE IS A DIFFERENT QUANTITY.** The meridional
streamfunction applied one constant density — `r_air`, the SURFACE value — at every level of a
mass-flux integral, i.e. it labelled a volume flux kg/s. Corrected, `Psi_max` goes 851.68 ->
373.84 (1e9 kg/s): **2.28x overweight at the cell core**, location unmoved. This is the tree
that uses that diagnostic to judge cell strength (`24ff23a` "revive Hadley/Ferrel cells",
`1e59daa` jet spin-down).

## Constants that are THIS TREE'S and must not be copied from ATHAD

| | here | ATHAD ships |
|---|---|---|
| `zeta` (radial stretch) | 3.715 | 3.0 |
| `checkRadialMetric` spread | **23.21x** — the worst in the family | 11.8x |
| `ATM_GRID_PTOP` | **0.08538** | 1e-6 |
| `ATM_GRID_BETA` | **3.988** | 4.33 |
| shell / scale height | **2.01** | 5.06 |

**`ATM_GRID_PTOP = 1e-6` on Earth's column sits at 81.9 km.** Importing ATHAD's default would
not regrid this model, it would **extend the shell 5.1x** while looking like a grid option. The
same goes for `buildReferenceColumn`: ATHAD integrates a DRY ADIABAT because its column is one
by construction, and that law puts this tree's 16 km lid at 117 K.

## The ported A/B knobs, and what each MEASURED here

All default to what this tree has always done, and all are verified bit-identical when unset.
**Every one behaved differently from ATHAD — porting its DEFAULTS would have been wrong four
times out of four.**

| knob | here | ATHAD |
|---|---|---|
| `ATM_PROJ_SWEEPS` | **inert**: -0.04 % at 10x, -0.07 % at 100x, cost 8x and 96x | -52.5 % |
| `ATM_METRIC_EXACT` | null on every INTEGRATED quantity — but see below | 2.8x worse |
| `ATM_BUOY_TREF` / `_CONSISTENT` | unmeasured, and the budget says not worth it | 5.49x at the surface |
| `ATM_GRID_PRESSURE` | ~1 %, and the sign is against it | +61 % on its free branch |
| `ATM_RAD_TOPO` | **NEW HERE, and it is this tree's defect, not ATHAD's** | inapplicable — no topography |
| `ATM_RHIE_CHOW` | ported, **unmeasured here** | -2.55x on the zonal Nyquist |

`ATM_PRESS_SWEEPS` exists too and is deliberately separate from `ATM_PROJ_SWEEPS`: this
projection calls the solver 200 times, so one knob would make "10 sweeps in the time loop" also
mean 2000 relaxations at startup, and any comparison would differ in its INITIAL STATE as well
as in the quantity under test. ATHAD lost an attribution exactly that way.

## Open risks

- **THE RADIATION COLUMN STARTED AT LEVEL 0, WHICH IS THE GROUND ONLY OVER OCEAN**
  (`ATM_RAD_TOPO`, default off, 2026-08-27). `MultiLayerRadiation` had
  `const int i_mount = 0;  // surface / bottom layer`. Over topography levels
  `0 .. i_topography-1` are rock, and `ThermoAtm`'s barometric loop writes a real `p_stat` into
  every one of them — so they carry MASS, `dp` is positive through the rock and enters `sum_dp`,
  and since each layer's optical depth is `tau_dry*dp_i/sum_dp`, **every air layer over every
  mountain was diluted by the rock beneath it**. `tau_layer`/`tau_above` were wrong ABOVE the
  ground, not merely inside it; the surface energy balance, `e_surf`, `epsilon_2D` and `T_air1`
  were all placed at sea level. **This is the family's Earth-constant pattern with TERRAIN in
  place of a number.** Off-branch bit-identical. With it on at 4 iterations: `epsilon_2D` max
  relocates from Angola to **29N 88E, the Himalaya** (0.088 -> 0.119), `tau_layer` max +25 % and
  also to the Himalaya, `tau_above` max -2.7 %, `radiation` min 95.2 -> 111.2 W/m2.
  **Recommended for the default after a longer run — not flipped, because the band constants
  here are tuned.**
- **`brunt_N2`'s +-0.03 to 0.047 s^-2 "boundary-layer" extrema were THE TERRAIN** — the centred
  difference at the first air level straddled it. Repaired unconditionally; ocean columns
  bit-identical. Extremes are now -0.00068/+0.0051, and the old ones sat on the Andes
  (14S 71W) and the Himalaya (29N 87E). **`ba5f542`'s note that they appear "every time" is
  retired.**
- **`p_dyn` CARRIES A GRID-SCALE MODE, AND HERE IT IS A SPIN-UP TRANSIENT** (2026-08-27). The
  Poisson operator is the compact 7-point Laplacian at `dr` while `div_src` and the Step-3
  gradient correction are `2*dr` central differences, which annihilate the Nyquist mode exactly.
  Measured here at pressure solve: per-axis Nyquist share **0.443 / 0.331 / 0.660**, absolute
  ~1.6e-06 on all three — **near-isotropic**, unlike ATHAD's k-only 0.961.
  **IT DOES NOT DECAY, AND THE RATIO SAID OTHERWISE.** Over 100 iterations the global index falls
  0.63 -> 0.0384, which was first written up here as a spin-up transient. **That was wrong.** The
  denominator, rms `p_dyn`, grows **24x** over the same run and the ABSOLUTE amplitude RISES:
  2.990e-06 (diagnostic 1) -> 3.735e-06 (13) -> **4.454e-06** (25). The per-axis measure agrees
  once the first transient is past — zonal Nyquist 4.34e-05 -> 1.98e-06 by diagnostic 5, then
  back UP to 3.62e-06, while the zonal anomaly grows faster, so the share falls 0.34 -> 0.19.
  **The stripes leave the plot because the smooth field grows around them, not because the noise
  dies.** *Caveat*: neither instrument isolates the Nyquist mode from smooth curvature — both are
  second differences — so a growing smooth field contaminates them; the share is the more
  trustworthy. The global index is normalised by rms `p_dyn` and **cannot be compared between
  trees**, because that denominator differs by seven orders; read the per-axis absolutes.
- **`Psi` DOES NOT CLOSE AT THE GROUND, IT IS WRITTEN BY THE INITIAL CONDITION, AND THE ATHAD
  LEVER IS ALREADY MEASURED INERT HERE** (2026-08-27, the first runs with `Psi` as a field).
  Measured by item 68's method — **`Psi(ground)` as an RMS over latitude, never as a max**:

  | iter | `Psi(ground)` rms | `Psi(lid)` rms | ground_rms / max\|Psi\| |
  |---|---|---|---|
  | 1 | 1.5629e+11 | **0.0000e+00** | 0.407 |
  | 3 | 1.5771e+11 | 0.0000e+00 | 0.414 |
  | 20 | 1.6225e+11 | 0.0000e+00 | 0.427 |
  | 100 | 1.6657e+11 | 0.0000e+00 | 0.401 |

  **The lid closes exactly and the ground does not.** Together those say the column-integrated
  meridional mass flux does not vanish. **It is NOT a boundary-cell artefact**: rms\|Psi\| decays
  smoothly through the whole depth — 1.666e+11 at 0 m, 1.598e+11 at 236 m, 1.269e+11 at 2163 m,
  7.05e+10 at 6088 m, 3.29e+07 at 14567 m — so it is a column-wide offset, not a spike at `i = 0`.

  **At iteration 1 it is already 93.8 % of its iteration-100 value** (ATHAD item 68 measured
  98.7 %), so it is written by the initial velocity profile and the dynamics only add 2.7 % over
  the next 99 iterations.

  **AND `Psi_max` MIGRATES ONTO THE DEFECT.** At iterations 1-3 it is at 45N / 1729 m — the real
  circulation — and by iteration 20 it is at 15N / **0 m** and stays there to 100. The interior
  circulation decays (this tree's jet spin-down) while the ground offset does not, so the same
  scalar reports the circulation early and the defect later. **Any `Psi_max` comparison that
  straddles that crossover is comparing two different quantities.**

  **THE INITIAL PROJECTION IS EXONERATED, AND IT IS `VelocityInitializer`** (`ATM_PSI_PROJ_DUMP=1`
  writes the streamfunction either side of `project_initial_velocity` as iterations -1 and -2):

  | stage | `Psi(ground)` rms | `Psi(lid)` rms |
  |---|---|---|
  | before the projection | 1.5543e+11 | 0.0000e+00 |
  | after the projection | 1.5541e+11 | 0.0000e+00 |

  **-0.0143 %.** The projection does not remove the non-closure; it does not touch it. So the
  offset is written by `VelocityInitializer`'s analytical Hadley/Ferrel profile and survives
  intact. **This also explains why `ATM_PROJ_SWEEPS` is inert here** where ATHAD got 52.5 % from
  it (-0.04 % at 10x, -0.07 % at 100x in the knob table above): upstream the projection cannot
  see this mode at all, so sweeping harder changes nothing.

  **AND IT CONTRADICTS THE COMMENT AT THE CALL SITE**, which says the projection "lets the time
  loop start from a clean div v = 0 state with the solenoidal part of the prescribed flow
  intact". For this quantity that is false.

  **AND `u` AT BOTH WALLS IS EXACTLY ZERO, WHICH CLOSES THE ARGUMENT.** `BC_Atm.h:679` sets
  `u.x[im-1] = 0` and `:697` sets `u.x[0] = 0` (with `un` in lockstep), both unconditional -- the
  latter carries a comment recording that leaving it unconstrained drove +-100 m/s coastal
  blowups. Verified in the written field at iteration 100, not merely in the assignment:
  `max|u|` = **0.000000e+00** at `i = 0` and at `i = im-1`, against ~7e-03 one cell in.

  So: for a zonally symmetric field with no mass flux through top or bottom, `div(rho v) = 0`
  forces the column-integrated meridional flux constant in latitude and zero at the poles, hence
  `Psi(ground) = 0`. Both walls are shut, and it is 1.55e+11 anyway. **The field is not
  divergence-free, and this is item 72 / item 86's non-adjoint projection with a physical
  consequence attached at last**: the residual divergence integrates to a **spurious net
  meridional mass transport, ~40 % of the Hadley cell's own strength**.

  **THE CAUSE IS THE PRESCRIBED PROFILE, AND `ATM_V_MASSBAL=1` REMOVES 94.8 % OF IT AT
  INITIALISATION AND 71.2 % AT ITERATION 100** (default off). `VelocityInitializer::init_v_or_w`
  builds `v` as a **linear ramp in height** from `coeff_sl` at the surface to `coeff_trop` at the
  tropopause, then a linear decay to the lid. Nothing constrains `INT(rho*v*dz) = 0`, or even
  `INT(v*dz) = 0`: for a linear ramp the volume integral is `H*(v_s+v_t)/2`, zero only if the two
  hand-set endpoints are exact opposites, and the mass integral needs something different again
  because rho decays roughly exponentially while the ramp is linear in z. At 15N the column runs
  +3.67 m/s at the ground against -0.54 m/s at 7.4 km — the poleward branch in the dense lower
  6 km, the return in thin air — and `INT(v*dz)` is **+9151 m^2/s**, so it fails to close in
  VOLUME as well. Two defects stacked, and the density weighting would still be wrong if the
  endpoints were opposites.

  **This is the family's rho-blindness one file upstream of where it was already caught**:
  `MinMax_Atm.cpp`'s `const double rho = r_air` was the DIAGNOSTIC version (`dabbc94` in ATHAD,
  `ede4810` here). The instrument was corrected; the thing it measures never was.

  The repair is the initial-condition analogue of ATHAD's `initBalancedState` — impose the
  constraint rather than hope the projection removes it. Per fluid column subtract the
  density-weighted column mean, `v <- v - INT(rho*v*dz)/INT(rho*dz)`, which zeroes the column
  mass flux while shifting the profile by a single constant, so the SHEAR defining the cell is
  untouched. 65 341 columns, largest correction 2.6e-01 non-dim.

  | iter | ground rms off | ground rms on | change | max\|Psi\| off | max\|Psi\| on |
  |---|---|---|---|---|---|
  | init | 1.5543e+11 | 8.0210e+09 | **-94.8 %** | | |
  | 20 | 1.6225e+11 | 3.1714e+10 | -80.5 % | 3.7984e+11 | 2.0371e+11 |
  | 100 | 1.6657e+11 | 4.7935e+10 | **-71.2 %** | 4.1568e+11 | 1.6077e+11 |

  **IT IS NOT A CURE AND THE RESIDUAL IS RISING MONOTONICALLY** — 8.02e+09 at initialisation to
  4.79e+10 at iteration 100, still climbing. **No limit is claimed; this file has been caught
  extrapolating that shape.** So the initial profile is the DOMINANT source but not the only one:
  the dynamics regenerate the offset, which is what a projection that converges to a
  non-divergence-free fixed point would do. Two further caveats: the residual is ~5 % even at
  initialisation, because `Psi` uses the ZONAL MEAN `rvbar` whose fluid-cell count varies with
  height where `i_topography` varies with `k`, so a per-column balance is not exactly a
  zonal-mean balance; and `Psi_max` still migrates onto the ground by iteration 100, though at
  iterations 20-80 it sits at 3.7-5.5 km — **the real circulation, which the control never
  showed after iteration 3.**

  **That makes `Psi(ground)` an INTEGRATED instrument for the projection residual** -- and a far
  more legible one than the local `div(u)` rms, which is what this tree and ATHAD have been
  reading. It is also why `ATM_PROJ_SWEEPS` cannot help: sweeping harder converges the solver to
  the same non-divergence-free fixed point item 72 measured as bit-identical under 64x. Operator weights differ too: **`c_phi/c_r` = 0.0322 here against
  0.59 in ATHAD**, a 16 km shell over a 6370 km radius against a 300 km one — so "the horizontal
  directions are weakly constrained" is plausibly true HERE and was measured false there.

- **The meridional streamfunction does not close at the ground, and it is not the solver.**
  `Psi(ground)` must be zero; its RMS over latitude is **42 % of the interior maximum**, and
  `ATM_PROJ_SWEEPS` at 10x and 100x moves it 0.04 % and 0.07 % with the knob verified connected.
  In ATHAD roughly half that non-closure was projection convergence and the rest structural;
  **here it looks like all-structural**, which makes this tree a CLEANER specimen of ATHAD item
  72's non-adjoint collocated stencil (Rhie-Chow named as the un-done repair) than ATHAD itself.
- **`ATM_METRIC_EXACT` halves the vertical wind and nothing integrated notices.** RMS ratios
  0.611/0.537 on two orthogonal slices, **p50 ratios 0.460/0.597**, every level 0.50-0.69 — the
  median falls as far as the max, so it is not a clipped extremum. `Psi` is built from the
  MERIDIONAL wind and KE from the horizontal components; a RADIAL metric governs the vertical
  one, and neither instrument reads it. **Unattributed** between the exact Jacobian (21.8x
  larger at the surface) and the curvature term added with it (`-curv*df/dr`, `curv = zeta`);
  the discriminating arm — exact Jacobian, curvature forced off, ~12 min — HAS NOT BEEN RUN.
  This is why the default is off, and the first write-up of it was WRONG and is corrected in
  the README.
- **No `div(rho u)/rho` diagnostic.** ATHAD prints it every iteration and judges the projection
  by it; here closure has to be computed from `meridional_streamfunction_*.csv`. Porting that
  print is the obvious next instrument.
- **Read `Psi(ground)` as an RMS over latitude, never as a max.** The max sits inside one cell,
  so a change elsewhere reads as bit-identical while the field moves. ATHAD item 68's trap.
- **Thread count still changes results in the last digits.** OpenMP reduction order is not
  associative; a tropical precipitation probe differs 4e-6 between 2 and 4 threads. State the
  thread count with any number quoted from here.

## The build hazard, because it produced a crash that looked like a success

`-MMD` tracks headers only for objects ALREADY compiled with it. `cli/atm.o` was four weeks
stale through an entire port, and it holds `cAtmosphereModel model;` as a **stack local**: adding
members to the class made `main` reserve the OLD `sizeof` while `libatom.a` constructed the NEW
one, so the constructor ran off the end of main's frame onto the stack canary. Every run aborted
with `*** stack smashing detected ***` AFTER printing its results, with no compiler warning, and
**AddressSanitizer found nothing** — an ODR size mismatch is not an out-of-bounds access to
anything it tracks. It was found by disassembling `main` (object at `rbp-0x1600`, canary at
`rbp-0x18`, object `0x1678` long).

Hence `$(LIB_OBJ) ... : Makefile` in the Makefile. Do not remove it.

**`make clean` deletes `python/pyatom.cpp`, which is generated but TRACKED** — restore it with
`git checkout` or regenerate with `make python`.

## Relationship to the family

Siblings: `ATHAD` (Hadean, 250 bar, water-vapour-dominated), `ATHAD_COND` (its
post-condensation epoch), `ATHAD_PERID` (its magma-ocean epoch), and the giants `ATJUP`,
`ATSAT`, `ATURAN`, `ATNEPT`, plus `ASTIM`. The giants share none of this tree's atmosphere code
— no `exp_rm`, no `MoistConvection`, no `m_node_weights`.

**Fixes flow both ways, and the forks find defects this tree cannot.** Everything in the
*Correctness fixes* README section was found in ATHAD, where Earth's constants stop being true.
Three defects came back the other way on 2026-08-26 — the OCEAN's copies of the Poisson and
`residuum_old` races, and `MinMax_Hyd`'s thread-order-dependent tied extrema — because ATHAD has
no ocean to have found them in.

**Before re-deriving anything, check the unmerged branches.** The density-weighted
streamfunction was fixed on `atom-metric-fixes` (2026-08-14), never merged, invisible, and
implemented a second time twelve days later. A cross-reference is not a check, and neither is a
branch nobody points at.
