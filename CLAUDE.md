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

All but `ATM_V_MASSBAL` default to what this tree has always done, and all are verified
bit-identical when unset. **`ATM_V_MASSBAL` BECAME THE DEFAULT ON 2026-08-28, SO EVERY OTHER ROW
IN THIS TABLE WAS MEASURED AGAINST AN INITIAL STATE THAT IS NO LONGER THE DEFAULT** — the arms
are being re-run against the new one; `ATM_V_MASSBAL=0` restores the old branch exactly.
**Every one behaved differently from ATHAD — porting its DEFAULTS would have been wrong four
times out of four.**

| knob | here | ATHAD |
|---|---|---|
| `ATM_PROJ_SWEEPS` | **inert**: -0.04 % at 10x, -0.07 % at 100x, cost 8x and 96x | -52.5 % |
| `ATM_METRIC_EXACT` | null on every INTEGRATED quantity — but see below | 2.8x worse |
| `ATM_BUOY_TREF` / `_CONSISTENT` | unmeasured, and the budget says not worth it | 5.49x at the surface |
| `ATM_GRID_PRESSURE` | ~1 %, and the sign is against it | +61 % on its free branch |
| `ATM_RAD_TOPO` | **NEW HERE, DEFAULT ON since 2026-08-28** — this tree's defect, not ATHAD's | inapplicable — no topography |
| `ATM_RHIE_CHOW` | **null on `Psi(ground)`, +0.005 %** — see below | -2.55x on the zonal Nyquist |
| `ATM_V_MASSBAL` | **NEW HERE, AND DEFAULT ON SINCE 2026-08-28**: -94.8 % of `Psi(ground)` at init, -71.2 % at iter 100 | not ported yet |
| `ATM_CONV_ADJ` | **NEW HERE**: surface lapse -19.45 -> -9.76 K/km | ATHAD's own file, default off there too |
| `ATM_ANELASTIC` | ported, **null on `Psi(ground)` (-0.006 %)** — and the reason is structural, below | on by default there |
| `ATM_BUOY_MOIST` | **NEW HERE**: +4.9 % on `ubud_buoy`, but ONLY with `ATM_BUOY_CONSISTENT` | not ported |
| `ATM_PROJECT_IN_LOOP` | **NEW HERE**: null on `Psi(ground)` at 10 AND at 200 sweeps (-0.013 %) | no equivalent |

**`ATM_PROJECT_IN_LOOP=<sweeps>` HAS NOW BEEN RUN TO 200 SWEEPS (2026-08-28), AND SWEEPS ARE NOT
THE LEVER.**
Default 0 = off. It is a real velocity projection inside the time loop: seed `aux` from `u/v/w`,
relax the Poisson, apply `v <- v - grad(p)`, with `p_dyn` saved and restored around the call
because the loop's own pressure is live and read by the next RK4 stage.

**It is NOT "call `run()` and subtract the gradient", and that distinction is the useful part.**
In the time loop `aux_*` does not hold a velocity: `RHS_Atm_Turb.cpp:1140` sets
`aux_u = rhs_u + dpdr_exp`, a momentum TENDENCY. So the loop's `run()` computes the pressure that
makes the ACCELERATION divergence-free — a legitimate fractional-step variant, since if `div(u)`
starts at zero and every tendency is divergence-free then `div(u)` stays zero — and subtracting
THAT gradient from `u` would be wrong. **So "the velocity is never projected" is true but the
model is not simply omitting a step; it projects the tendency instead.**

Four arms at iteration 100, 24 threads:

| sweeps | `Psi(ground)` rms | vs off | `div(u)` at solve 25 | wall |
|---|---|---|---|---|
| off | 1.66571090e+11 | — | 1.766e-03 | 527 s |
| 10 | 1.66569245e+11 | -0.0011 % | 1.764e-03 | 556 s |
| 10 + anelastic | 1.66561990e+11 | -0.0055 % | — | 552 s |
| **200** | **1.66549029e+11** | **-0.0132 %** | **1.761e-03** | **1017 s** |

**THE 10-SWEEP NULL WAS DISCOUNTED FOR A REASON THAT NO LONGER HOLDS.** The first write-up said
the projection was "barely acting" because `project_initial_velocity` runs **200** passes and
this was given **10**, and red-black Jacobi converges the LARGEST scales slowest — error decaying
like `(1 - c/N^2)` per sweep, so a domain-scale mode on a 181x361 grid wants O(N^2), and
`Psi(ground)` is exactly such a mode. It has now had 200. **Twenty times the sweeps bought 2.5x
the `div(u)` reduction** (-0.28 % against -0.11 %) and left `Psi(ground)` where it was. That is
the same shape as `ATM_PROJ_SWEEPS` being bit-identical under 64x: **the solver converges to a
fixed point that is not divergence-free, and more relaxation converges to it harder.**

**THE COST ESTIMATE IN THE FIRST WRITE-UP WAS ALSO WRONG.** "~20x the solver cost per iteration"
is true per solve, but the solve is about a fifth of the runtime: 1017 s against 527 s, **1.93x
wall clock, 17 minutes.** The arm was deferred as expensive and was not.

**WHAT IT STILL DOES NOT SETTLE, AND THE MISSING INSTRUMENT IS ONE PRINT.** `div(u)` is only ever
reported AT THE PRESSURE SOLVE, never immediately after `project_velocity_in_loop`. So
"the projection zeroes `div(u)` and one time step puts all of it back" and "the projection never
reduces `div(u)` at all" have the IDENTICAL signature in this table, and they are different
defects. Do not read this null as either one.

**And the whole line may be chasing a quantity the model does not constrain**: `INT(v dz)` is
7.83e3 m^2/s, so `div(u)` was never zero to begin with, and `Psi(ground)` may be measuring
accumulated solver error rather than one defect. The volume-vs-mass explanation remains
**UNSUPPORTED** — neither confirmed nor refuted by anything run so far.

**`ATM_BUOY_CONSISTENT` IS NO LONGER UNMEASURED, AND THE ROW ABOVE SAYING "the budget says not
worth it" WAS WRONG ON ITS SECOND HALF** (2026-08-27). Measured: `ubud_buoy` **0.000003 ->
1.208688**, a factor of 400 000, which puts it **ninety-fold above `ubud_pgf`** (0.0135). ATHAD's
extra-`dt` defect (its items 34/42) is confirmed here at the same order of magnitude. What is
still true is that it does not move `Psi`: that force is RADIAL and `Psi` is meridional.

**`ATM_BUOY_MOIST`, and why it needs the row above.** The shipped buoyancy is
`(t - t_ref_level[i])`, TEMPERATURE ONLY. Water vapour is lighter than dry air, and the model
already computes the virtual temperature — `r_humid = p/((1 + R_W_R_A_m1*c - cloud - ice)*T)` —
uses it for density and the streamfunction, and the buoyancy discards it, because the Boussinesq
form takes a TEMPERATURE anomaly where a DENSITY anomaly belongs. The knob adds
`T_v = T(1 + 0.608 q_v - q_c - q_i)` at all four sites, with the 0.608 taken as
`R_WaterVapour/R_Air - 1` so it matches `r_humid` rather than being a second constant.
**The REFERENCE is virtual too** (`tv_ref_level`, built in the same sweep): `t_ref_level`'s whole
purpose is that the body force has zero mean at every height, so moistening the parcel alone
would have added a uniform updraft rather than a buoyancy.
**On the shipped branch it is a null in the 10th digit — because the buoyancy itself is inert.**
With `ATM_BUOY_CONSISTENT=1` it is **+4.9 %** on `ubud_buoy`. *Adding moisture to a force that is
not acting cannot show anything*, and that is the whole content of the first measurement.

**AND SURFACE EVAPORATION THEREFORE CANNOT DRIVE CONVECTION.** It moistens levels 0-3 and reaches
the momentum equation through nothing on the shipped branch. It also **never writes `t`**, so it
removes no latent heat and cannot cool the surface into stability either — the same shape as
`Q_Sensible`, which is written and read by nothing. Both halves of the surface energy exchange
are open, in opposite directions.

**THE VOLUME-VS-MASS EXPLANATION IS NOW MEASURED: REAL IN DIRECTION, FAR TOO SMALL TO BE THE
CAUSE** (2026-08-28). The arm that was said to be unreachable was reachable all along -- the TIME
LOOP never applies the pressure to the velocity, but `project_initial_velocity` DOES, and
`ATM_PSI_PROJ_DUMP=1` writes `Psi` either side of it. Measured across the initial projection with
`ATM_V_MASSBAL` on, so all arms start from `Psi(ground)` = 8.0210e+09:

| projection enforces | relaxations | after | change |
|---|---|---|---|
| `div(u)` — volume | 200 | 8.0038e+09 | -0.215 % |
| `div(rho u)` — MASS | 200 | 7.9860e+09 | **-0.436 %** |
| `div(rho u)` — MASS | 4 000 | 7.8258e+09 | -2.434 % |
| `div(rho u)` — MASS | 12 000 | 7.7834e+09 | -2.962 % |
| `div(rho u)` — MASS | 40 000 | 7.7671e+09 | **-3.166 %** |
| `div(u)` — volume | 40 000 | 7.8873e+09 | **-1.667 %** |

**MASS BEATS VOLUME BY 1.9x AT CONVERGENCE** (2.03x at 200 relaxations, 1.90x at 40 000), so the
explanation in `PressureSolverAtm.h:129` is confirmed in DIRECTION. **AND IT ACCOUNTS FOR 1.5
PERCENTAGE POINTS OF A 100 % NON-CLOSURE.** Even the mass projection, fully converged, leaves
96.8 % of `Psi(ground)` standing. Volume-vs-mass is a real contributor and NOT the cause.

**THE PROJECTION IS CONVERGED, AND THE MIDDLE OF THAT CURVE LIES.** Increments are +2.00, +0.53,
+0.20 points for 200 -> 4 000 -> 12 000 -> 40 000. Read at 4 000 alone it looks convergence-limited
and this file briefly said so; the plateau near -3.2 % is unambiguous by 40 000. **So item 72's
reading stands: the solver converges to a fixed point that is not divergence-free**, and ~97 % of
the non-closure is structural, outside any projection's reach. 40 000 relaxations is O(N^2) on a
181x361 grid, so this is not a sweep-count question.

**WHY `ATM_ANELASTIC` IS A NULL ON A FULL RUN, AND WHY THAT IS NOT A RESULT ABOUT CONTINUITY.** The port is
faithful — source `div(u*) + u*_r dln(rho_bar)/dr`, the matching operator term folded into the
existing `num_a` off-diagonal, base state built in `densities()` at line 514 and so available to
`project_initial_velocity` at line 554. It moves `Psi(ground)` **-0.006 %**. The reason is the
structural gap recorded above: in the time loop `PressureSolverAtm::run()` computes `p_dyn` and
**never applies `v <- v - grad(p)`**, so changing which continuity the PRESSURE solves for cannot
change the VELOCITY during a run. **Consequence for the volume-vs-mass diagnosis: it is
UNSUPPORTED, not confirmed.** The two-integrals measurement stands on its own, but the arm that
would have tested the explanation could not reach the velocity. Testing it needs the velocity
projected INSIDE the time loop — a change to how the model steps, not a knob.

`ATM_RHIE_CHOW` being a null on `Psi(ground)` is STRUCTURAL, not a failure to tune: `D4`
annihilates smooth fields by construction — that is the property that makes the knob safe — and
`Psi(ground)` is a domain-scale quantity. **So the collocated-grid checkerboard and the
streamfunction non-closure are NOT the same defect**, despite both having been attributed to
div/grad non-adjointness. At smooth scales the 2*dr and compact operators agree to O(dx^2), which
cannot produce a 40 % error.

`ATM_PRESS_SWEEPS` exists too and is deliberately separate from `ATM_PROJ_SWEEPS`: this
projection calls the solver 200 times, so one knob would make "10 sweeps in the time loop" also
mean 2000 relaxations at startup, and any comparison would differ in its INITIAL STATE as well
as in the quantity under test. ATHAD lost an attribution exactly that way.

## Open risks

- **THIS TREE HAD NO DRY CONVECTIVE ADJUSTMENT, SO NOTHING COULD REMOVE A SUPERADIABATIC LAYER**
  (`ATM_CONV_ADJ`, default off, ported from ATHAD 2026-08-27). No such file, no call site, only
  `MoistConvection.h`. Measured over 100 iterations at 24 threads: at iteration 0 the surface
  layer is STABLE (-3.93 K/km, `brunt_N2` < 0 in **0.000** of columns); by iteration 20 the
  lowest 39 m runs at **-18.05 K/km, twice the dry adiabat**, with `brunt_N2` < 0 in 59 % of
  columns at i = 0 and i = 1 and 0 % at i = 2 — one unstable layer, exactly the
  surface-to-first-level step. It is manufactured by the TIME LOOP: the surface warms +0.76 K
  over the run while the air above warms +0.15 K, and nothing mixes them.
  **With the adjustment on**: surface lapse **-19.45 -> -9.76 K/km** against a -9.8 dry adiabat,
  `brunt_N2` median at i = 0 **-2.69e-04 -> -1.0e-07** and at i = 1 -4.75e-05 -> **+6.09e-05**
  (frac<0 0.587 -> **0.006**), i = 3 and above unchanged to four figures. **Enthalpy drift
  3.07e-16**, machine precision; 67.3 % of columns adjust, worst column 2 sweeps of an allowed
  64. `Psi(ground)` -0.02 %, `max|Psi|` -0.3 %, radiation extremes +0.11 % / -0.55 %.
  **READ THE MAGNITUDE, NOT THE FRACTION**: `frac<0` at i = 0 stays 0.590 while the value falls
  to -1e-07 — that residual is the `tol_nd` tolerance on a layer that is now neutral by
  construction, so the fraction counts round-off there. **And max/min temperature are
  BIT-IDENTICAL between the arms** (35.065051 / -37.203979) while the median surface profile
  moves 0.3 K — the extremum is a false null, again.
- **`ATM_SFC_COUPLED` WAS WRITTEN, MEASURED AND REMOVED, and the lesson is worth more than the
  knob.** It solved the surface energy balance together with a neutral-lapse constraint, to cure
  the same `brunt_N2` < 0. It never fired: at setup the column is stable, so the guard correctly
  declined, and the two arms were identical at iteration 0. **`MultiLayerRadiation` is called
  seven times, ALL at lines 713-889, while the iteration loop starts at line 985 — MLR never runs
  inside the loop and cannot maintain anything at iteration 100.** Check where a routine is
  CALLED before attributing a standing feature to it.
  **THAT PARAGRAPH IS WRONG, AND IT BROKE ITS OWN RULE** (2026-08-28). Lines 711-892 are LAMBDA
  DEFINITIONS — `refresh_radiative_teq`, `apply_radiative_heating`, `refresh_radiation_diag`,
  `cloud_radiation_diag` — and the CALL SITES are at **1463-1471, inside the loop**. The active
  `radiation_mode` is **5** (the run log prints it), which calls `cloud_radiation_diag()`, and
  that runs MLR, **every `teq_refresh_stride` = 20 iterations**. So MLR does run in the loop, on a
  20-iteration cadence — which is exactly the cadence of the `epsilon` step measured between
  iterations 20 and 40. The sentence above told the next reader to check where a routine is
  CALLED, and was itself written from where one is DEFINED.
  **`ATM_SFC_COUPLED`'s "it never fired" needs re-examining on that basis**; the guard-declined-
  at-setup half of it stands, the "MLR cannot maintain anything" half does not.
- **`brunt_N2 < 0` AT `i = 0` IS THE OCEAN MASK, AND IT IS A BOUNDARY-CONDITION GAP**
  (2026-08-28). Measured at iteration 100 on the level-0 radial slice:

  | | cells | `P(N2 < 0)` | median `N2` |
  |---|---|---|---|
  | ocean | 43 602 (66.73 %) | **0.9998** | **-3.70e-04** |
  | land | 21 739 (33.27 %) | **0.0049** | **+2.11e-04** |

  43 593 of the 43 700 unstable cells are ocean and 107 are land. **The "59 % of columns"
  recorded above was always the ocean fraction** — the number was in the file for a day and
  nobody matched it to the mask.

  Four facts, and the last is the defect:
  1. **Level 0 is not prognostic**: `RungeKutta_Atm_Turb.cpp:111` is `for(int i = 1; i < im-1;)`.
  2. **`t` has NO radial BC at `i = 0`.** `bcRadius` carries three pattern lists —
     `both_cubic`, `vn_bot_cubic_top`, `cubic_bot_vn_top` — and `t` is in **none** of them,
     while `p_stat`, `r_humid`, `cloud`, `ice` and `tke` all get one. `t` has explicit special
     handling at the LID (`pin_t_top`) and nothing at the bottom.
  3. **The mechanism meant to cover that is a self-assignment over ocean.** `BC_Atm.h:433`,
     documented as "copy surface-level (`i_mount`) values down to the i=0 reference layer … so
     surface fluxes see the correct surface state", is `t.x[0] = t.x[i_topography]`. Over land
     `i_topography > 0` and it does real work — which is why **land is stable**. Over ocean
     `i_topography == 0` (`FileIO_Atm.cpp:488`) and it reads `t.x[0] = t.x[0]`.
  4. **So level 0's only coupling to the air above is a Shapiro smoother**, `damp_wiggles(t,…)`
     at `cAtmosphereModel.cpp:1113`, whose `along_i` pass starts at `i_surf`. Not a physical
     mixing process. Level 0 drifts **+0.786 K** over 100 iterations while the air above moves
     +0.15 K.

  **The stencil explains why it is ONE layer.** `N2` at `i = 0` is one-sided on `theta[0],
  theta[1]`; at `i = 1` it is centred on `theta[0], theta[2]`; at `i = 2` it is centred on
  `theta[1], theta[3]` — the first stencil that does not see level 0. **The 59/59/0 pattern is
  one bad value read by two stencils, not two unstable layers.**

  **THE MODEL HAS NOT DECIDED WHETHER OCEAN LEVEL 0 IS THE OCEAN OR THE AIR.** If it is a
  prescribed SST skin then `N2 < 0` there is PHYSICALLY CORRECT — a heated surface is
  superadiabatic against the air — and the repair is to apply the bulk flux to level 1. If it is
  the lowest air level it needs RK4 and a flux BC. It is currently neither.
  **This reframes `ATM_CONV_ADJ`**: its -19.45 -> -9.76 K/km is a palliative mixing away a
  boundary-condition gap, not a missing convection scheme.

  **THE ATTRIBUTION HAS BEEN RUN AND THE PARAGRAPHS ABOVE ARE WRONG ON THE MECHANISM**
  (`ATM_T0_ATTRIB=1`, 2026-08-29, default off and read-only: it snapshots levels 0 and 1 and
  differences them after every stage that can write `t`, cos-lat means in kelvin, cumulative,
  ocean and land separately). `nm = 100`, 24 threads, the default arm:

  | stage | ocean i=0 | ocean i=1 | land i=0 | land i=1 |
  |---|---|---|---|---|
  | `BC_Atm` | -0.0014 | -0.0940 | -30.8878 | +126.5407 |
  | `RungeKutta` | **0.0000** | **-16.7068** | 0.0000 | -206.2882 |
  | `teq_relaxation` | **+0.8605** | **+16.9771** | +31.5685 | +78.5385 |
  | every other stage | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
  | **NET** | **+0.8591** | **+0.1763** | +0.6807 | -1.2090 |

  **`apply_teq_relaxation` IS THE ONLY THING THAT WRITES OCEAN LEVEL 0.** It runs EVERY
  iteration in `radiation_mode` 5 (`cAtmosphereModel.cpp:1521`), loops from `i = 0`, and relaxes
  at `omega_teq` = **0.20 per iteration** — an e-folding of ~4.5 iterations. So "level 0 is
  unconstrained" is FALSE: it is Newtonian-relaxed onto `t_eq` harder than anything else in the
  model. `BC_Atm` contributes **-0.0014 K in 100 iterations** over ocean, which is the
  self-assignment measured rather than inferred, and `RungeKutta` is **exactly 0.0000**, which
  is `RungeKutta_Atm_Turb.cpp:111` measured rather than inferred. **`unattributed` is 0.0000, so
  the twelve hooks account for the whole field** — the instrument is closed, not sampled.

  **AND IT IS NOT A DRIFT. IT IS A ONE-TIME APPROACH THAT ARRIVES BY ITERATION 40.** Ocean i=0
  NET by checkpoint: **0.8153 / 0.8604 / 0.8588 / 0.8590 / 0.8591** at iterations 20 / 40 / 60 /
  80 / 100 — flat to four decimals for the last sixty iterations. The recorded "+0.786 K over
  100 iterations", read here and in the README as a secular drift, is a STEP: level 0 starts
  ~0.86 K below `t_eq.x[0]` and relaxes onto it. What put it below is the CO2 perturbation —
  `[co2-perturb DIAG]` prints the surface `t_eq` shift as **+1.075 K at build and 1.169-1.179
  thereafter**, i.e. the target moves +0.10 K over the whole run and level 0 tracks it. Ocean
  i=1 is the one that genuinely drifts, and only just: 0.1726 -> 0.1763, +0.0037 K per 20
  iterations, still rising.

  **THE SUPERADIABATIC LAYER IS `t_eq` RELAXED ONTO A LEVEL THE DYNAMICS CANNOT REACH.** Level 0
  has no RK4 term at all, so it equilibrates EXACTLY onto its target and stays there. Level 1
  has a dynamics term of **-16.71 K cumulative** that the relaxation (+16.98) must fight, so it
  settles BELOW its own target — two terms ~95x the +0.176 net. The step between the two levels
  is precisely the amount by which the dynamics hold level 1 off a target level 0 is free to sit
  on. **So the defect is not a missing boundary condition and a Dirichlet pin would not remove
  it**: pinning level 0 to `t_eq.x[0]` is what the relaxation already achieves. The defect is
  that a surface-trapped CO2 perturbation is applied to `t_eq` at a level with no dynamics, and
  the level above it cannot follow.

  **TWO CORRECTIONS TO THE PARAGRAPHS ABOVE, AND THE SECOND IS THE FAMILY'S OWN RULE AGAIN.**
  (1) "Level 0's only coupling to the air above is a Shapiro smoother, `damp_wiggles(t)`" is
  wrong twice: the actual writer is the relaxation, and **`damp_wiggles(t)` DOES NOT RUN AT ALL
  in these runs** — it sits behind `moist_phys_active` and `moist_phys_start_iter` is **300**
  against `nm` = 100, so every measurement in this bullet, and every `brunt_N2` and
  `ATM_CONV_ADJ` figure recorded beside it, comes from a DRY run in which that smoother never
  fired. (2) The mechanism was written from the control flow rather than from a measurement, for
  the third time in this file. The instrument cost one build and two runs.

  **WHAT SURVIVES.** The ocean/land split, the stencil argument for why it is one layer, and the
  reframing of `ATM_CONV_ADJ` as a palliative all stand — it palliates a target inconsistency
  rather than a boundary-condition gap.

  **`ATM_SFC_FLUX` IS NOW MEASURED, IT IS A NULL, AND THE REASON IS THAT A 100-ITERATION RUN IS
  TWENTY SECONDS LONG.** `ATM_SFC_FLUX=15`, `nm = 100`, 24 threads, against the arm above:
  ocean i=1 `RungeKutta` **-16.7068 -> -16.6464** (+0.36 %), `teq_relaxation` 16.9771 -> 16.9209,
  NET **0.1763 -> 0.1808**, i.e. **+0.0045 K on a 0.68 K step**. Ocean i=0 unmoved (0.8591 ->
  0.8590). The relaxation absorbs **93 %** of everything the flux adds, in the same iteration.

  **THE TIMESCALES, NOW PRINTED BY THE MODEL AT EVERY RUN** (`[TIMESCALES]`, unconditional —
  ATHAD's item 47 names the missing duration as an open risk and this tree inherited the gap):

  | | value |
  |---|---|
  | advective unit `metricShellLength()/u_0` | 16024 m / 8 m/s = **2002.9 s** |
  | one iteration, `dt_visc` = 1e-4 | **0.2003 s** |
  | `nm = 100` | **20.03 s of physical time** |
  | `nm = 400` (the shipped default) | 80 s |
  | `omega_teq` = 0.20 **per ITERATION** | e-folding **1.001 s** |
  | `ATM_SFC_FLUX` c_H = 15 W/m2/K | **3293 s (0.92 h)** |

  **`omega_teq` CARRIES NO `dt` AND THEREFORE NO TIMESCALE.** Every physical term in `rhs_t` is
  multiplied by `dt` and scales with the step; the relaxation is a fixed fraction per iteration
  and does not. At 0.20 it is a **1-second** process competing with a **3293-second** surface
  flux — **3290x faster, by construction rather than by tuning.** The run reaches
  `1 - exp(-20.03/3293)` = **0.605 %** of the skin-air difference; the measured NET change is
  0.0045/0.70 = **0.64 %**. Prediction and measurement agree to 6 %, which is what says the knob
  is correctly implemented and correctly scaled and simply has no time to act.

  **SO ITEMS 1 AND 2 ARE ONE DEFECT AND NEITHER IS REACHABLE BY A SURFACE KNOB AT THIS `nm`.**
  Nor by a larger one: one hour of physical time is 17 973 iterations and one day is 431 000, at
  ~20 s of wall clock each. **The lever is `dt`, or the relaxation, not the iteration count.**

  **AND IT EXPLAINS `radiation_mode` 5's OWN DESIGN.** The code comment at
  `cAtmosphereModel.cpp:733` records that the in-RHS Held-Suarez term at its PHYSICAL strength
  (`k_a` = 1/4 day) "enters as ~1e-8/iter, never imposes `t_eq`" — which is exactly this
  arithmetic — so it was replaced by a per-iteration numerical fraction to make the CO2 signal
  persist. **The relaxation is a workaround for the run being twenty seconds long, and it is now
  the thing that dominates every surface diagnostic in this tree.** Read `ATM_CONV_ADJ`'s
  -19.45 -> -9.76 K/km in that light: a convective adjustment mixes to neutral INSTANTANEOUSLY,
  with no `dt` and no timescale, which is why it is the only surface treatment that has ever
  moved this layer.

- **`Q_Sensible` IS WRITTEN AND READ BY NOTHING** (`RHS_Atm_Turb.cpp:513`), **AND IT IS NOT THE
  SURFACE FLUX — an earlier version of this bullet ran two different quantities together.**
  The array is `coeff_S * lap(T)`, the CONDUCTIVE heat-flux divergence, computed at every
  interior cell; adding it anywhere would double-count the thermal diffusion `diffusion_t_re`
  already applies. The surface bulk flux is a SEPARATE term, `c_H*(T_s - T_air1)` with
  `c_H = 15 W/m2/K` at `MultiLayerRadiation.h:425`, which the surface balance debits and nothing
  credits to the air. Both are real defects; only the second one bears on the instability above,
  and **MLR's copy of it is refreshed only every 20 iterations** (`radiation_mode` 5 ->
  `cloud_radiation_diag()` at `cAtmosphereModel.cpp:1471`; an earlier version of this bullet said
  MLR never runs in the loop at all, which was wrong) — so the flux still has to be formed in the
  loop from the live `t.x[0]`/`t.x[1]` if it is to act every step.

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
  **FLIPPED ON BY DEFAULT 2026-08-28, AND THE LONGER RUN IS WHY — THE 4-ITERATION COMPARISON WAS
  TAKEN BEFORE THE DEFECT APPEARS.** Off branch, `nm = 100`, level-0 radial slice:

  | iter | max `epsilon_2D` | at | max `tau_layer` | at |
  |---|---|---|---|---|
  | 20 | 0.08805 | 17S 12E — **Angola** | 0.12070 | 20S 68W |
  | **40** | **0.66300** | **28N 88E — Himalaya** | **7.59045** | **28N 88E** |
  | 100 | 0.67683 | 28N 88E | 7.88527 | 28N 88E |

  Between iterations 20 and 40 there is a STEP — `epsilon_2D` x7.5, `tau_layer` **x63** — and the
  maximum relocates to the Himalaya permanently. The recorded off-branch figure (0.088, Angola)
  reproduces exactly at iteration 20, so the old measurement was not wrong; **it was taken at 4
  iterations, before the thing the knob repairs exists.** At iteration 100 the off branch is
  **7.7x** the number the knob was judged against.

  **THE MECHANISM, READ OFF THE HIMALAYA COLUMN AT 28N (`k = 88`), ITERATION 100:**

  | lvl | `p_stat` | `Epsilon` | `tau_layer` | `q_v` |
  |---|---|---|---|---|
  | 0 | 951.2 | **0.67683** | **7.88527** | **7.295** |
  | 1 | 946.7 | 0.01778 | 7.88527 | **0.000** |
  | 2 | 941.7 | 0.01862 | 7.88527 | 0.000 |

  ocean at the same latitude: `Epsilon` 0.0217, `tau_layer` 0.0220, `q_v` 12.8.
  `tau_layer` is CONSTANT through the rock because of the fill at `MultiLayerRadiation.h:485`
  (applied on both branches), so the plotted value is the ground value replicated downward and
  the excess is in the ground value itself — 360x the ocean's. `Epsilon` is NOT filled off-branch,
  so its 0.677 is MLR's own level-0 value, 38x the cell above it.
  **The cause is `q_v` = 7.295 at level 0 and exactly 0.000 in every rock level above.** `BC_Atm`
  Pass 3 copies the mountain-top humidity down (`c.x[0] = c.x[i_mount]`) while the rock is dry, so
  in `tau_i = tau_dry*dp_i/Sum(dp) + tau_wv*vp_i/Sum(vp)` with `vp_i = c_i*dp_i`, level 0 holds
  essentially the whole sub-surface vapour path and collects a large share of the water-vapour
  optical depth.
  **`ATM_RAD_TOPO=0` restores the sea-level column exactly.** The band constants are still tuned
  on the old branch — that caveat stands and is now the argument for re-checking them, not for
  leaving the radiation at sea level over every mountain.
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

  **`Psi(lid)` = 0 IS THE INTEGRATION CONSTANT, NOT A RESULT, and the "lid closes / ground does
  not" framing gave it weight it never had.** `MinMax_Atm.cpp:232` integrates DOWNWARD from the
  lid with `psi[im-1] = 0` by construction, so that column of the table is arithmetic and is
  identical in every arm ever run. **All of the content is in `Psi(ground)`**, which is therefore
  identically the whole-column integral `2*pi*a*cos(phi)*INT(rho*vbar dz)` — not a value AT
  `i = 0` but a TOTAL, reported at `i = 0`. It does not vanish, and it must.
  **It is NOT a boundary-cell artefact**: rms\|Psi\| decays
  smoothly through the whole depth — 1.666e+11 at 0 m, 1.598e+11 at 236 m, 1.269e+11 at 2163 m,
  7.05e+10 at 6088 m, 3.29e+07 at 14567 m — so it is a column-wide offset, not a spike at `i = 0`.

  **AND IT CANNOT BE RELOCATED OR RESOLVED AWAY. TWO TEMPTING NON-CURES, BOTH RULED OUT**
  (2026-08-28):

  1. **Pushing the offset up the column by reshaping the initial `v(z)` cannot work.**
     `Psi(ground)` is a definite integral, so it is INVARIANT to any vertical redistribution of
     `v` that preserves `INT(rho*v dz)`. The profile can be reshaped so `Psi` is near-zero
     through most of the depth with the whole discrepancy in one thin layer, and the ground
     value will not move by one digit. There is no "where" for a total to live. The only initial
     condition that moves it is one that changes the column integral — which is `ATM_V_MASSBAL`.
     Re-integrating UPWARD from the ground instead is the cosmetic version: it relabels which end
     carries the constant, changes no velocity, and destroys the instrument.
  2. **Increasing `im` cannot close it either — the defect is in the INTEGRAND, not the
     quadrature.** Checked directly, by re-integrating the written field with a cubic spline on
     the SAME 41 levels (nodal `rvbar` recovered from the `psi` differences; the back-out is
     marginally stable, so the lid anchor was swept over its plausible range and the three
     anchors agree to three decimals):

     | field | trapezoid (model) | cubic spline | change |
     |---|---|---|---|
     | shipped, iter 20 | 1.622498e+11 | 1.619763e+11 | **-0.169 %** |
     | shipped, iter 100 | 1.665711e+11 | 1.663176e+11 | **-0.152 %** |
     | `ATM_PROJECT_IN_LOOP=200`, iter 100 | 1.665490e+11 | 1.662955e+11 | -0.152 % |
     | `ATM_V_MASSBAL=1`, iter 100 | 4.793468e+10 | 4.777011e+10 | -0.343 % |

     **Quadrature is ~0.15 % of a 40 % non-closure.** In ABSOLUTE terms it is ~1.6e8-2.5e8 kg/s
     in every arm — so even the residual `ATM_V_MASSBAL` leaves is **200x larger than the
     discretisation error**, and refining the grid cannot reach that either.
     The reason is visible in the source without running anything: `VelocityInitializer.h:55`
     sets the 15N Hadley column with `init_v_or_w(m.v, 75, -3.0, 4.0)` — tropopause **-3.0**,
     surface **+4.0**. For the linear ramp the CONTINUOUS integral is `H*(v_s+v_t)/2`,
     proportional to `+4.0 - 3.0 = +1.0`, nonzero before any grid exists. More levels converge
     the sum to that nonzero limit. **`ATM_V_MASSBAL` removing 94.8 % by changing `v` on the SAME
     grid with the SAME quadrature is the other half of the same proof.**
     *Caveat worth keeping*: this tree's stretch is the family's worst (`zeta` 3.715, spread
     23.21x, 38.9 m bottom layer), so a real trapezoid error could hide near the surface — the
     measurement above is what bounds it at 0.15 %, and it is a bound, not an argument.

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

  **THE CAUSE IS THE PRESCRIBED PROFILE, AND `ATM_V_MASSBAL` REMOVES 94.8 % OF IT AT
  INITIALISATION AND 71.2 % AT ITERATION 100** (**DEFAULT ON since 2026-08-28**; `=0` restores
  the old branch, verified at iteration 100: `Psi(ground)` 1.6657e+11 and `max|Psi|` 4.1568e+11,
  the recorded values to five figures). `VelocityInitializer::init_v_or_w`
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
- **`ATM_METRIC_EXACT` does not damp the vertical wind, it produces a DIFFERENT one from
  iteration 0.** RMS ratios 0.611/0.537 on two orthogonal slices, **p50 ratios 0.460/0.597**,
  every level 0.50-0.69 — which the first two write-ups read as a halving. It is a relocated
  structure, not a damped one. `Psi` is built from the MERIDIONAL wind and KE from the
  horizontal components; a RADIAL metric governs the vertical one, and neither instrument
  reads it.
  **THE ATTRIBUTION ARM HAS BEEN RUN, AND THE QUESTION HAS NO ANSWER** (`ATM_METRIC_NOCURV`,
  `6d4c68c` 2026-08-26 — an earlier version of this bullet said it had NOT been run, and that
  line was two days stale). Difference RMS against legacy, zonal / longitudinal: Jacobian only
  **92.7 % / 83.1 %**, curvature only **72.8 % / 65.3 %**, **both together 60.8 % / 51.8 %** —
  the combined change is SMALLER than either piece alone, so the two halves partially cancel
  and there is no additive split to attribute. The test's own premise failed too: the arms
  already differ **at iteration 0** by 95 % of the field, because `project_initial_velocity`
  uses `exp_rm` in the Poisson operator AND in the gradient correction, so the metric acts
  through the INITIAL PROJECTION first and no arm ever started from a shared state.
  **NOT DECIDABLE BY MORE A/B RUNS** — it is a physics judgement about which field is right,
  and the prerequisite is the `div(rho u)/rho` print in the next bullet, not another arm.
  This is why the default is off; the full record is in the README.
- **No `div(rho u)/rho` diagnostic.** ATHAD prints it every iteration and judges the projection
  by it; here closure has to be computed from `meridional_streamfunction_*.csv`. Porting that
  print is the obvious next instrument.
- **THE CLOUD OPTICAL DEPTH SATURATES ALOFT ON THE SHIPPED BRANCH, AND `ATM_CLOUD_TAU_MAX` IS
  DEFAULT ON SINCE 2026-08-28.** `cwp_cap_col` caps the COLUMN condensate path; it does not bound
  a LAYER, and the two are different constraints. Full-field `max epsilon` at iteration 100:

  | arm | max `epsilon` (all levels, all latitudes) | where |
  |---|---|---|
  | shipped, no ceiling | **0.999624** | 28N 88E, **4988 m** |
  | shipped + ceiling (**the default**) | **0.885490** | 20N 53E, 3316 m |
  | `ATM_RAD_TOPO=1`, no ceiling | 0.999747 | 28N 88E, 0 m |
  | `ATM_RAD_TOPO=1` + ceiling | **0.891525** | 35N 84E, 0 m |

  **No cell anywhere exceeds 0.892 with the ceiling on, on either branch.** A near-blackbody layer
  is the pathology the de-saturation split (`MultiLayerRadiation.h:169`) exists to prevent — it
  pins the emission level and collapses the OLR — and it was present in the SHIPPED model at
  ~5 km, with `ATM_RAD_TOPO` off. The topo branch does not create it; it moves it to the ground
  where a level-0 diagnostic finally showed it. Flipped on top-of-atmosphere evidence: clear-sky
  OLR **bit-identical** at 180.33882 W/m2, cloudy +0.053. `ATM_CLOUD_TAU_MAX=0` disables.
  The knob scales `LWP_i` AND `IWP_i` together so the LW `tau_cloud` and the SW albedo bump see
  the same condensate — capping one and not the other is the imbalance `cwp_cap_col`'s own note
  warns about.

- **`Results_Atm.cpp:40` HAS BEEN PRINTING THE FULL-FIELD `max epsilon`, WITH ITS HEIGHT, IN EVERY
  RUN LOG.** The table above took three rounds of slice extraction — level 0, then one longal cut,
  then a proposal to run 70N/70S arms — to reach a number already in the logs being grepped.
  **Third occurrence of this pattern in the family**: ATHAD item 68 ("every number in this item
  came from re-reading VTK files that had been sitting in `output_Hadean/` unexamined") and item
  72 ("the model has been printing this residual in every run log all along"). **A level-0
  diagnostic cannot answer a column question**, and reaching for the existing print costs nothing.

- **THE "OLR" WAS THE LID TEMPERATURE, NOT AN OUTGOING FLUX — AND THE GREENHOUSE WAS NEVER 60 W/m2
  WRONG** (2026-08-28). `cloud_radiation_diag` printed `radiation.x[im-1]`, the TOP LAYER's own
  emission, as "cos-lat-mean top radiation". It read ~180 W/m2 against Earth's ~240 and was cited
  all day as a 60 W/m2 model error. It is **sigma*T_lid^4**: the lid sits at **236.15 K** and
  `sigma*236.15^4 = 176.3 W/m2`. The instrument was measuring the lid pin.

  **The proof is one line of the new print**: zeroing cloud takes the column optical depth from
  **28.10 to 2.18** — twenty-six optical depths — and moves that "OLR" by **2.7 W/m2**. A real OLR
  moves by tens.

  `column_olr()` integrates the upward flux properly, as the mirror of the `L_down` sum MLR
  already does — each layer's emission attenuated by the layers above, plus the ground through the
  whole column, starting at `i_topography` so it never integrates through rock:

  | | model | Earth |
  |---|---|---|
  | column tau, clear | **2.18** | ~1.5-2 |
  | column tau, cloudy | **28.10** | — |
  | **OLR clear** | **273.89** | ~265 |
  | **OLR all-sky** | **193.31** | ~240 |
  | **cloud LW forcing** | **80.59** | **~25** |

  **THE CLEAR-SKY RADIATION IS ESSENTIALLY RIGHT**, 273.9 against ~265, which confirms the
  de-saturation split's own claim of "OLR ~263" and vindicates `eps_dry`, the Bignami 0.0056 and
  `co2_band_scale`. **THE WHOLE DEFECT IS CLOUD**: 26 of the 28 optical depths, and a longwave
  forcing **3.2x** Earth's. `cwp_cap_col = 250 g/m2` is still 2.5x the observed ~100, and there is
  **no cloud fraction at all** — every column is treated as fully overcast.

  **AND THIS REVERSES TWO THINGS SAID EARLIER THE SAME DAY.** "Re-tuning `cwp_cap_col` would not
  move the OLR" was based on the broken instrument and is **wrong** — it is now the single lever
  that matters. And the standing objection to `ATM_RAD_TOPO`, that it costs 0.54 W/m2 of OLR in a
  model already 60 low, is **void**: that 0.54 was a change in the lid layer's emission, and the
  model is not 60 low. **Fourth instrument-shaped defect in this tree** after `Psi`'s constant
  density, `Q_Sensible`, and `brunt_N2`'s terrain extrema: computed, printed, trusted, and
  measuring something other than its name.

- **THE CLEAR-SKY OLR EXCESS IS THE TEMPERATURE PROFILE, NOT THE BAND CONSTANTS — AND THE TWO
  ERRORS HAVE OPPOSITE SIGNS** (2026-08-28). Clear-sky OLR is 273.9 against Earth's ~265, a
  +9 W/m2 excess. It decomposes into two larger errors that partly cancel:

  | | OLR clear | vs Earth |
  |---|---|---|
  | constants on a CORRECT column (`test/rad_selftest`, US-standard, CO2 380) | **255.0** | **-10** |
  | the model's own atmosphere | **273.9** | **+9** |
  | difference — the profile | **+18.9** | |

  **So the constants UNDER-produce on a correct column** — 255.0 against Earth's ~265 and against
  the de-saturation note's own claim of "~263" — and the model's profile more than cancels it.
  Tuning `eps_dry` or `co2_band_scale` to remove the +9 would push a column that is already
  10 W/m2 too transparent further wrong.

  **THE PROFILE BIAS GROWS WITH HEIGHT AND PEAKS AT THE EMISSION LEVEL**, model minus US-standard
  at 28N, iteration 100:

  | z | 0 m | 612 m | 2160 m | 6079 m | **9908 m** | 13220 m |
  |---|---|---|---|---|---|---|
  | model - std | +7.67 | +7.59 | +11.43 | +15.76 | **+18.75** | +19.57 |

  Lapse rate **5.38 K/km against 6.50** standard. The photosphere sits at ~9.9 km (`tau_above` = 1),
  which is exactly where the bias is largest, so the column emits from a level **18.75 K too warm**.
  **AND THE RADIATION SCHEME ITSELF PRODUCES THAT BIAS**: the offline harness, given a US-standard
  column, returns **241.41 K at 9923 m against the 223.7 K it was handed** — +17.7 K, matching the
  full model's +18.75 almost exactly. It is the scheme's own radiative equilibrium, not something
  the dynamics did to it. **`ATM_CONV_ADJ` WAS THE OBVIOUS LEAD AND IS NOW A MEASURED NULL ON IT**
  (`nm = 100`, same binary, current default): tropospheric lapse **5.381 -> 5.362 K/km**, 0.35 %,
  and the profile above 82 m is unchanged to **+0.002 K at every level** — including +0.002 K at
  9908 m, the emission level carrying the whole +18.75 K bias. Clear-sky OLR moves
  273.899 -> 273.791, **-0.11 W/m2**. The adjustment touches level 0 only (-0.191 K).
  **The reason is structural and was predictable**: the scheme mixes where `brunt_N2 < 0`, i.e.
  SUPERadiabatic layers, and this troposphere has the opposite defect — at 5.38 K/km it is far
  MORE stable than the 6.5 moist or 9.8 dry adiabat, so there is nothing for it to act on. A
  convective adjustment cannot steepen a profile that is already too stable.
  **So the warm-aloft bias is the radiation scheme's own equilibrium**, which the offline harness
  shows independently: handed a US-standard column it returns +17.7 K at 9923 m unaided.

  **AND IT IS NOT REACHABLE BY THE CLEAR-SKY CONSTANTS. SWEPT IN THE HARNESS**
  (`ATM_EPS_DRY`, `ATM_CO2_BAND`, both default to the shipped values, ~0.07 s per case):

  | case | T_sfc | T(9.9 km) | lapse K/km | OLR |
  |---|---|---|---|---|
  | **US-standard / Earth** | **288.15** | **223.75** | **6.49** | **~265** |
  | shipped 0.684 / 0.17 | 280.76 | 241.41 | 3.97 | 255.0 |
  | `eps_dry` 0.50 | 276.46 | 242.49 | 3.42 | 271.0 |
  | `eps_dry` 0.95 | 292.66 | 234.50 | 5.86 | 182.5 |
  | `co2_band` 0.0 | 261.50 | 243.15 | 1.85 | 251.4 |
  | `co2_band` 1.40 | 293.23 | 232.52 | 6.12 | 165.2 |

  **LAPSE RATE AND OLR ARE LOCKED IN OPPOSITION, AND EARTH IS OFF THE CURVE.** Every setting that
  steepens the lapse toward 6.5 collapses the OLR: 6.12 K/km costs **165 W/m2**. Earth has BOTH
  6.49 and 265, and no combination of these two constants comes near it. The constants move the
  profile's OFFSET — the surface spans 261-293 K — and barely its SLOPE aloft, where 9.9 km stays
  in 232-243 K against a standard 223.75 across the entire sweep. **Tuning them is exhausted.**

  **THE STRUCTURAL READING, and it closes the `ATM_CONV_ADJ` null above.** This scheme computes
  RADIATIVE equilibrium; Earth's troposphere is in radiative-CONVECTIVE equilibrium. Earth's
  radiative equilibrium is strongly UNSTABLE near the ground (~15 K/km), which is what triggers
  convection and lets it set 6.5. This scheme's is **3.97 K/km — far too STABLE**, so there is
  nothing for a convective adjustment to mix, which is exactly what `ATM_CONV_ADJ=1` measured.
  **The lead is the optical-depth weighting**: `tau_i = tau_dry*dp_i/Sum(dp)` distributes by MASS
  alone, with no pressure broadening, whereas real LW absorption goes roughly as `p*dp` and so
  concentrates optical depth near the surface — which is what makes a real radiative equilibrium
  steep at the bottom. Mass-only weighting spreads it too evenly and flattens the profile.
  **TESTED AND REFUTED** (`ATM_TAU_PBROAD`, exponent on `(p_i/p_0)`, default 0 = shipped and
  bit-identical, column total renormalised so it redistributes rather than adds):

  | pbroad | T_sfc | T(9.9 km) | lapse K/km | OLR |
  |---|---|---|---|---|
  | 0 (shipped) | 280.76 | 241.41 | 3.97 | 255.0 |
  | 0.50 | 282.30 | 241.97 | 4.06 | 271.0 |
  | 1.00 | 283.45 | 242.87 | 4.09 | 282.7 |
  | 2.00 | 284.98 | 244.30 | 4.10 | 296.3 |

  The lapse moves **3.97 -> 4.10 K/km across the whole range, 3 %**, and the upper troposphere gets
  **WARMER** (241.4 -> 244.3), the opposite of the prediction. The knob is plainly connected — OLR
  swings 255 -> 296 — it simply does not act on the profile's slope.

  **AND THAT PAIRS INTO SOMETHING SHARPER: THE TOTAL OPTICAL DEPTH MOVES THE LAPSE, ITS VERTICAL
  DISTRIBUTION DOES NOT.** `co2_band` 0 -> 1.4 swings the lapse 1.85 -> 6.12; `pbroad` 0 -> 2 pins
  it at 4.0-4.1. That is backwards for a radiative-transfer scheme, where the profile follows
  `tau(z)` directly. Checked against the grey analytic form `sigma*T^4 = A(1 + 3*tau/2)`,
  normalised at the lid:

  | level | `tau_above` | `sigma T^4` | grey prediction | ratio |
  |---|---|---|---|---|
  | 40 | 0.000 | 176.3 | 176.3 | 1.000 |
  | 35 | 0.814 | 196.1 | 391.7 | **0.501** |
  | 30 | 3.851 | 277.1 | 1195.0 | **0.232** |
  | 20 | 22.148 | 376.9 | 6034.4 | **0.062** |

  **The temperature barely responds to the optical depth it sits under.** *Caveat, and it matters*:
  that table is the FULL MODEL's `t`, which dynamics also shape, so a departure from pure
  radiative equilibrium is partly expected. What is not explained by dynamics is the OFFLINE
  harness — no dynamics at all — returning the same flat 3.97 K/km from a US-standard input.
  **THE TRIDIAGONAL SOLVE HAS NOW BEEN READ AND TESTED, AND IT IS EXONERATED ABOVE 3.7 km AND
  GUILTY BELOW IT — AND THE SCHEME HAS NO RADIATIVE EQUILIBRIUM AT ALL** (`test/rad_iterate`,
  2026-08-29, ~1 s per pass, no model change).

  **1. THERE IS NO FIXED POINT.** `rad_selftest` calls the scheme ONCE and reports what comes
  back; that answers "does it run", not "what is its equilibrium". Applied REPEATEDLY to its own
  output on a US-standard column it DIVERGES:

  | pass | T(0 m) | T(2163) | T(6088) | T(9923) | lapse 0-9.9 km | OLR |
  |---|---|---|---|---|---|---|
  | input | 288.15 | 274.09 | 248.58 | 223.65 | 6.500 | — |
  | 1 | 280.76 | 280.05 | 260.04 | 241.41 | **3.965** | 255.01 |
  | 2 | 268.31 | 281.23 | 264.31 | 247.76 | 2.070 | 256.57 |
  | 3 | 243.81 | 278.49 | 264.37 | 249.31 | **-0.554** | 240.87 |
  | 4 | 193.67 | 270.63 | 259.89 | 246.35 | -5.309 | 207.82 |
  | 5 | **NaN** | 254.36 | 249.10 | 237.64 | NaN | NaN |

  The surface collapses at an ACCELERATING rate (-7.4, -12.5, -24.5, -50.1 K) and the lapse
  inverts by pass 3. **So "the radiation scheme's own radiative equilibrium" names something that
  does not exist**, and the 3.97 K/km is ONE STEP of a divergent map applied to a US-standard
  column — not a converged state. Every constant ever fitted against a single pass was fitted
  against a transient. (Not yet separated: how much of the divergence is the tridiagonal and how
  much the surface-balance override bolted on beside it.)

  **2. AGAINST ANALYTIC GREY RE ON THE SCHEME'S OWN `tau_above`, THE SOLVE IS RIGHT ALOFT.**
  `sigma*T^4 = (F/2)(1 + 3*tau/2)` with F = the scheme's own OLR — not fitted, not normalised, so
  the lid value is a PREDICTION:

  | z [m] | 16023 | 9923 | 6088 | 3678 | 2163 | 1211 | 613 | 82 | 0 |
  |---|---|---|---|---|---|---|---|---|---|
  | `tau_above` | 0.000 | 0.355 | 0.686 | 1.014 | 1.305 | 1.547 | 1.739 | 1.961 | 2.003 |
  | model | 219.52 | 241.41 | 260.04 | 272.43 | 280.05 | 283.92 | 284.45 | 277.40 | 280.76 |
  | grey | 217.76 | 242.29 | 259.91 | 274.38 | 285.57 | 293.95 | 300.14 | 306.83 | 308.06 |
  | **model - grey** | +1.76 | -0.88 | +0.13 | -1.95 | **-5.52** | **-10.04** | **-15.68** | **-29.43** | -27.30 |

  Within ~2 K over the whole column above 3.7 km, including an unfitted lid. **The tridiagonal
  solve is doing grey radiative transfer correctly where the optical depth is thin.**

  **3. BELOW 3.7 km IT FAILS, AND IN THE LOWEST 800 m IT INVERTS THE PROFILE.** The model peaks
  at **284.64 K at 819 m** and falls DOWNWARD to **277.40 K at 82 m**, where grey RE requires a
  monotone rise to 308.06 K at the ground. That is 16 K of error OF THE WRONG SIGN in the bottom
  800 m, and the surface-balance override then puts 280.76 K underneath it — **a superadiabatic
  surface step manufactured by the radiation scheme itself**, in a harness with no dynamics.

  **4. SO TWO CONCLUSIONS RECORDED ABOVE ARE WRONG.**
  **(a) The +18.75 K warm-aloft bias is NOT a solver defect.** Grey RE on the scheme's own tau is
  **+18.64 K above US-standard at 9.9 km** (242.29 against 223.65) and the scheme returns
  241.41 — it reproduces grey RE there to 0.88 K. The bias is the ABSENCE OF CONVECTION and
  nothing else, which is what a radiative-equilibrium scheme is supposed to give.
  **(b) "This scheme's radiative equilibrium is 3.97 K/km — far too STABLE, so there is nothing
  for a convective adjustment to mix" IS AN ARTEFACT OF THE BOTTOM FAILURE.** The scheme's OWN
  grey equilibrium over 0-2163 m is (308.06 - 285.57)/2.163 = **10.40 K/km — SUPERADIABATIC**,
  against a 9.8 dry adiabat. There is a great deal for a convective adjustment to mix; it is
  hidden by a solver that flattens the bottom 3.7 km. **`ATM_CONV_ADJ`'s measured null was taken
  on a profile made flat by a solver defect, not by physics, and the convective route to Earth's
  6.5 K/km is REOPENED.**

  **THE WHOLE MATRIX IS WRONG, NOT ITS BOTTOM, AND POINT 4(a) ABOVE IS RETRACTED**
  (`ATM_RAD_EQUIL`, default 0, 2026-08-29). The grey comparison in point 2 used the Eddington
  **3/2 diffusivity factor** and this scheme has none, so it was the wrong yardstick, and it
  exonerated the solve aloft when it should not have. Two REFERENCE-FREE tests replace it.

  **(i) THE OUTPUT DOES NOT SATISFY THE BALANCE THE SCHEME IS WRITTEN TO SOLVE.** For an air
  layer, absorbed = emitted is `eps_i*(U_i + D_i) = 2*eps_i*B_i`, i.e. **`B_i = (U_i + D_i)/2`
  with `eps` cancelling** — no optical-depth assumption, no diffusivity factor, nothing external.
  Sweeping U and D from the scheme's OWN epsilon and OWN output:

  | z [m] | 16023 | 9923 | 6088 | 2163 | 1368 | 613 | 82 |
  |---|---|---|---|---|---|---|---|
  | residual W/m2 | +3.76 | -17.69 | -45.52 | -76.48 | **-79.97** | -72.69 | -28.06 |
  | as dT [K] | +1.55 | -5.75 | -12.25 | -16.81 | **-16.95** | -15.08 | -5.99 |

  The output is too WARM against its own balance at every level below 10 km, peaking at
  **-80 W/m2 (-17 K) at 1.4 km**. Which is also why iterating it DIVERGES: the balance wants it
  colder, and each pass overshoots.

  **(ii) SOLVED TO CONVERGENCE, THE SAME BALANCE GIVES A COMPLETELY DIFFERENT PROFILE.** Jacobi
  on the scheme's own epsilon and shortwave, 200 000 sweeps, validated by a TOA closure of
  **+0.000 W/m2** that is enforced nowhere in the iteration:

  | z [m] | 16023 | 9923 | 6088 | 2163 | 819 | 82 | 0 |
  |---|---|---|---|---|---|---|---|
  | scheme | 219.52 | 241.41 | 260.04 | 280.05 | 284.64 | 277.40 | 280.76 |
  | **correct** | **193.91** | **207.03** | **218.47** | **235.72** | **244.22** | **250.48** | **269.85** |
  | diff | -25.61 | -34.38 | -41.57 | -44.32 | -40.41 | -26.93 | -10.91 |

  Wrong by **25-44 K at EVERY level**, not only below 3.7 km. And the correct answer aloft is
  **16.6 K COLDER than US-standard** at 9.9 km, not warmer — which is what radiative equilibrium
  should be, since convection warms the upper troposphere relative to it. **So the "+18.75 K
  warm-aloft bias" IS a solver defect after all**, and 4(a)'s exoneration was an artefact of the
  3/2 yardstick.

  **(iii) POINT 4(b) STANDS AND IS STRENGTHENED.** Correctly solved, the lower troposphere is
  **15.78 K/km over the bottom 2 km** from the shipped shortwave and **19.93 K/km** from a
  realistic one — hugely superadiabatic against a 9.8 dry adiabat. This file's own structural
  reading ("Earth's radiative equilibrium is strongly UNSTABLE near the ground (~15 K/km), which
  is what triggers convection... This scheme's is 3.97 K/km — far too STABLE") had the physics
  exactly right and the attribution exactly wrong: **it is the solver, not the scheme's physics.**
  `ATM_CONV_ADJ`'s null was measured on a profile the solver flattened.

  **THE REPLACEMENT IS EXACT, AND CHEAPER THAN WHAT IT REPLACES.** `U[i]` crosses the interface
  BELOW layer i and `D[i]` the one ABOVE it, so `U[i] - D[i]` is not a net flux at one level —
  the net at layer i's upper interface is `U[i+1] - D[i]`, and radiative equilibrium holds THAT
  constant at F. Eliminating B and D from the transfer step gives a single downward sweep:

      U[i+1] = U[i] - eps_i*F/(2 - eps_i)
      B_i    = U[i] - F/2 - eps_i*F/(2*(2 - eps_i))

  from `U` above the lid = F = absorbed shortwave, with `B` at the ground = `U[i_mount+1]`
  because a blackbody surface does not attenuate its own emission. The `eps/(2 - eps)` is the
  layer's opacity to its own emission; assuming `eps/2` costs 2.23 K, which is how it was caught.
  **Verified against the 200 000-sweep Jacobi at max |dT| = 0.000000 K.** O(im) and
  non-iterative, against the shipped O(im^2) `CC` construction.
  **And it gives the scheme a fixed point**: `ATM_RAD_EQUIL=1` converges and is EXACTLY
  stationary thereafter (`max|dT|` = 0.0000), where the shipped path reaches NaN in five passes.
  Off-branch bit-identical — `rad_selftest` reproduces 280.76 / 241.41 / 255.0.

  **IT MUST NOT BE FLIPPED ON, BECAUSE IT UNCOVERS A SECOND DEFECT: THE SHORTWAVE IS ~2.5x TOO
  LOW AND THE BROKEN SOLVER WAS MASKING IT.** `rad_equator_short` = **163.3 W/m2** against
  Earth's annual-mean equatorial TOA insolation of ~416 (global mean 341), `rad_pole_short` =
  100.0 against ~173. The shipped scheme emits 255 W/m2 while absorbing 150 — **a 105 W/m2 TOA
  imbalance that was invisible because the solve was never energy-closed.** Correct the solver
  alone and the ice-albedo feedback runs away:

  | equatorial SW | converges | surface | lapse 0-2.2 km | lapse 0-9.9 km | OLR |
  |---|---|---|---|---|---|
  | 163.3 (shipped) | pass 3 | **219.12 — snowball** | 12.81 | 5.14 | 65.3 |
  | **416 (Earth)** | **pass 1** | **340.92** | **19.93** | **7.998** | **382.7** |

  340.9 K with a strongly superadiabatic lower troposphere is the textbook Manabe-Strickler
  radiative equilibrium. **Two errors have been cancelling: a solver biased warm and a shortwave
  forcing 2.25x too weak.** Fixing either alone makes the model worse, which is the family's
  `mue_ch4` pattern (ATNEPT item 2) a third time.

  **THE SHORTWAVE IS NOW FIXED TOO, AND THE PARABOLA WAS THE WRONG SHAPE RATHER THAN THE WRONG
  SCALE** (`ATM_SW_INSOL=<S0>`, default 0 = off, bit-identical; 1361 = Earth). `parabola(x) =
  x^2 - 2x` on `j/j_half` is a quadratic in LATITUDE; the annual-mean insolation is exactly a
  quadratic in SIN(latitude),
  `S(phi) = (S0/4)*(1 - 0.477*P2(sin phi))`, where the 0.477 is the obliquity's second Legendre
  coefficient — **solar geometry, not a fit**. Cos-lat means:

  | | equator | pole | global mean | rms vs true |
  |---|---|---|---|---|
  | true annual-mean | 421.4 | 178.0 | **340.3** | — |
  | shipped 163.3 / 100.0 | 163.3 | 100.0 | **151.3** | 198.3 |
  | true endpoints in the parabola | 416 | 173 | 370.0 | 39.7 |
  | best-fit parabola | 404.9 | **63.6** | 340.3 | 16.7 |

  **The two constants cannot both be right in the parabolic form**: the true endpoints overshoot
  the integral by 8.7 % because the parabola is too flat in mid-latitudes, and a least-squares
  parabola gets the integral right only by driving the pole to 63.6 W/m2 against a true 178. In
  the correct form there are no free constants at all, only `S0`. `rad_equator_short` /
  `rad_pole_short` are ignored on this branch and restore the old profile exactly when it is off.

  **THE FOUR ARMS, AND WHY THE TWO KNOBS ARE A PAIR** (US-standard column, iterating the scheme on
  its own output, `test/rad_iterate`):

  | arm | pass 1 | fixed point | lapse 0-9.9 km | TOA closure |
  |---|---|---|---|---|
  | shipped | 280.76 | **NaN by pass 6** | 3.965 | never closes |
  | `ATM_SW_INSOL` only | 297.08 | **140.92, still falling** | 5.490 | — |
  | `ATM_RAD_EQUIL` only | 269.85 | 219.12 — snowball | 6.331 | +0.000 |
  | **BOTH** | **342.02** | **342.02, stationary from pass 2** | **8.024** | **+0.000** |

  **FIXING THE SHORTWAVE ALONE IS WORSE THAN SHIPPING NEITHER** — it diverges to 140.92 K with an
  INVERTED lapse of -11.7 K/km. Together they converge in ONE pass, are exactly stationary
  thereafter, close the top of atmosphere to +0.000 W/m2 with closure enforced nowhere, and give
  a lower troposphere of **19.998 K/km over the bottom 2 km**. That is clear-sky
  Manabe-Strickler radiative equilibrium, and it is what a convective adjustment is supposed to
  find and mix. **`MultiLayerRadiation::warnIfHalfRepaired()` prints a loud warning if exactly
  one of the two is set**, because that arm is the worst of the four and looks like a repair.

  **STILL OPEN, AND NOT CHASED HERE**: the scheme has no atmospheric shortwave absorption and no
  Rayleigh scattering, so with `albedo_equator` = 0.1 it absorbs 387.7 W/m2 at the equator
  against Earth's ~316 — the albedo constants are a SURFACE albedo doing a PLANETARY albedo's
  job. That is the next thing in this chain, and it is why 342.02 K should be read as "the
  scheme's clear-sky RE, correctly solved" rather than as a claim about Earth.
  *Side note, not the same quantity*: the SURFACE lapse at 28N is **194.65 -> 188.74 K/km** — a
  7.6 K jump across the 39 m from level 0 to level 1, the boundary-condition defect above, which
  the adjustment dents by 3 %. Do not confuse it with the tropospheric lapse; they move
  independently and only the tropospheric one bears on the OLR.

- **`test/rad_selftest` SEGFAULTED FROM 2026-08-26 TO 2026-08-28 AND NOBODY RAN IT.**
  `MultiLayerRadiation` gained `tau_above` / `tau_layer` with the ported instruments and also
  reads `i_topography` and `short_wave_radiation`; none of the four were allocated in the harness,
  so the first write past an unallocated `Array` took the process down (exit 139). **The offline
  reference column — the thing the band constants were fitted against — was unrunnable for the
  whole period in which those constants were being argued about**, including a full afternoon of
  this session. Repaired by allocating everything the scheme touches.

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
