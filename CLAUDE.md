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
| `ATM_RAD_TOPO` | **NEW HERE, AND STILL DEFAULT OFF** — this tree's defect, not ATHAD's. The two "flipped on" claims below are WRONG, corrected 2026-09-02: `MultiLayerRadiation.h:293` is `e && atoi(e) != 0`, and every run's own `[RUN CONFIG]` banner prints `RAD_TOPO=0*` | inapplicable — no topography |
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

## The ice scheme was ThreeCat for one day, and is TwoCat again since 2026-09-01

**`CategoryIceScheme` = 2 (TwoCat) IS THE DEFAULT AGAIN**, in `param.py` and in all four tracked
configs, at the user's instruction on 2026-09-01 and for the reason in the subsection at the end
of this section: **ThreeCat's precipitation is a clamp residual, not a sum of rates.** Set 3 to
get ThreeCat back; its three repairs are written or named and none is finished. The rest of this
section is the ThreeCat record, kept because it is what the revert was decided on.


`CategoryIceScheme` = **3** in `param.py` and in all four tracked configs (was 2). The user had
asked for this some days earlier and it had not been carried out; `param.py` still had the
ThreeCat line commented out above a live TwoCat line.

**"ThreeCat is BROKEN (NaN)" IS STALE.** Run from the accepted configuration's own iteration-600
checkpoint, 600 -> 700 at 24 threads: exit 0, **zero NaN**, `max|w|` **17.97 m/s** against the
+-100 clamp, 10 min 27 s. `TwoCatIceScheme.h:189` still carries a comment citing "the ThreeCat
NaN blow-up" as the reason for a floor it added; whatever that was, it does not reproduce here.

**IT OVER-PRODUCES BY 3.8x AND THE EXCESS IS ALMOST ENTIRELY SNOW**, from the same checkpoint at
iteration 700:

| | TwoCat | ThreeCat | NASA |
|---|---|---|---|
| **Precipitation** | ~1200 | **3749 mm/a** | **978** |
| P_rain | 738 | 1056 | |
| **P_snow** | **0.29** | **2470** | |
| P_graupel | 0 | 464 | |
| P_conv | 0.23 | 0.049 | |
| max cell, all species | 1.8e+04 | **9.46e+04 mm/a** (259 mm/d) | |

**So the flip is a 3.8x regression on the model's headline output and it was made deliberately,
at the user's instruction, with this measurement in front of it.** The snow is the whole story:
2470 against 0.29 mm/a. `P_graupel` alone (464-920 mm/a depending on the print parity) is half
of NASA's total precipitation. Every species pegs the same 9.46e+04 mm/a per-cell ceiling, so a
cap is binding somewhere in all three.
**IDENTIFIED 2026-09-01**: it is `IceSchemeCommon::P_max_flux` = 3.0e-3 kg/(m2 s), and it is
not merely binding -- it is what sets the answer. See below.
**THE SNOW FIX IS NOW PORTED (2026-08-31), AND IT CUTS THE SNOW BY 62 %.** `c4d1bd1`'s two
changes — `c_c_au` 4.0e-4 -> 1.0e-3 and the SNOW-side riming `c_rim` -> `c_rim/5` — applied to
`ThreeCatIceScheme`. Warm-side shedding keeps the full `c_rim`, exactly as in TwoCat. Six
iterations from the same checkpoint:

| mm/a | before the port | ported | NASA |
|---|---|---|---|
| P_rain | 1190 | 1118 | |
| **P_snow** | **2751** | **1047** (-62 %) | |
| P_graupel | 921 | 919 (untouched) | |
| **total** | **4862** | **3084** (-37 %) | **978** |

**AND THE SAME TREATMENT APPLIED TO THE GRAUPEL RIMING** (`S_g_rim`, which TwoCat has no
equivalent of), six iterations from the same checkpoint:

| mm/a | before any port | + snow port | **+ graupel port** | NASA |
|---|---|---|---|---|
| P_rain | 1190 | 1118 | **1118** | |
| P_snow | 2751 | 1047 | **1047** | |
| P_graupel | 921 | 919 | **578** (-37 %) | |
| **total** | **4862** | **3084** | **2743** | **978** |

Rain and snow are unchanged to four figures, so the term is cleanly separated. Graupel falls
only 37 % rather than 80 % because riming is one of several graupel sources — `S_g_agg`,
`S_g_shed` and deposition are untouched.

**STILL 2.8x NASA, AND THE SNOW FRACTION IS THE NEXT THING.** Rain 1118, snow 1047, graupel 578:
snow is **38 % of the total** against Earth's 5-10 %, where the same fix took TwoCat to 16 %. So
the port helped and did not finish the job in ThreeCat.

**THE COLD-SIDE GRAUPEL RIMING WAS ALSO USING THE SNOW CONSTANT, AND THAT IS NOW FIXED — BUT
GRAUPEL IS NOT RIMING-LIMITED.** `c_g_rim` = 4.43 exists, is named for graupel riming, and was
used ONLY by the warm-side shedding `S_g_shed`; the cold-side `S_g_rim` reached for `c_rim` =
18.6, **4.2x larger**. The two corrections were measured SEPARATELY so neither is confounded:

| arm, 6 iterations from iteration 600 | P_rain | P_snow | P_graupel | total | x NASA | snow frac |
|---|---|---|---|---|---|---|
| before any port | 1190 | 2751 | 921 | 4862 | 4.97 | 56.6 % |
| + snow port | 1118 | 1047 | 919 | 3084 | 3.15 | 33.9 % |
| + graupel `c_rim/5` | 1118 | 1047 | 578 | 2743 | 2.80 | 38.2 % |
| **+ graupel `c_g_rim/5`** | **1118** | **1046** | **515** | **2679** | **2.74** | **39.0 %** |
| *Earth / NASA* | *~880* | *~70* | *~0* | *978* | *1.00* | *5-10 %* |

**A FURTHER 4.2x CUT IN THE RIMING COEFFICIENT BUYS 11 %** (578 -> 515). The constant is now the
right one and the change is free, but it settles the question: **the remaining graupel is not
riming, it is `S_g_agg`, `S_g_shed` and deposition**, none of which the port touches. That was
already visible in the first graupel arm falling 37 % where an 80 % fall was the naive
expectation.

**WHAT IS LEFT IS SNOW.** 1046 mm/a, **39 % of the total** against Earth's 5-10 %, and the
riming fix that took TwoCat to 16 % has already been applied. ThreeCat has snow sources TwoCat
lacks — `S_s_agg`, and the graupel/snow exchanges — and none of them has been looked at.
Precipitation is 2.74x NASA and the composition is wrong in the same direction: rain 42 %,
snow 39 %, graupel 19 %, against roughly 90 / 7 / 0.

**AND ONE MORE GAP, NOT CLOSED: `ThreeCatIceScheme`'s AUTOCONVERSION IS STILL GRID-MEAN.** It
reads `S_c_au = c_c_au*(q_c - 0.0002)` — a hardcoded 0.2 g/kg threshold on the GRID MEAN, where
TwoCat/OneCat/ZeroCat were all made fraction-aware in `3ea78e3` and now use `ATM_QC_CRIT`
(default 0.05 g/kg) against the IN-CLOUD value `q_c/f`. Under the now-default `ATM_CLOUD_FRAC`
the grid-mean `q_c` peaks near 0.06 g/kg, so `q_c - 0.0002` is negative almost everywhere and
ThreeCat's autoconversion is effectively OFF — its rain comes from accretion and melting. This
is the fifth occurrence of the grid-mean defect and the one module that never got the fix.
**CLOSED 2026-09-01** -- ported, and it is worth -598 mm/a; see the subsection below.

**`CategoryIceScheme` = 3 in the config selects ThreeCat**, and every measurement in this file
dated on or before 2026-08-31, or after 2026-09-01, was taken on TwoCat. The ThreeCat arms are
the ones dated 2026-08-31 evening and 2026-09-01, and they are labelled.

### The snow question was the wrong question: EVERY ThreeCat number is a clamp residual

**THE SNOW FRACTION WAS TO BE CHASED THROUGH `S_s_agg` AND THE GRAUPEL/SNOW EXCHANGES, ONE
COEFFICIENT PER RUN. THE INSTRUMENT ANSWERED IT IN ONE RUN AND THE ANSWER IS THAT NO COEFFICIENT
IN THIS SCHEME IS THE LEVER** (2026-09-01). `ATM_SS_DIAG=1`, new, print-only, default off: the
per-term budget of all three fluxes, walked exactly as the scheme integrates them, with the
`max(0, ...)` floor, the `P_max_flux` cap and the temperature window charged to their own bucket
so that `sources - sinks + clamp = P_x(ground)` is an IDENTITY. It closes to 1e-3 mm/a against
numbers of 1e+10, and to 1e-12 once the fluxes are physical.

Six iterations from the accepted configuration's iteration-600 checkpoint, 24 threads, cos-lat
global means in mm/a, on the SHIPPED scheme:

| snow | `S_i_au` | `S_s_agg` | `S_s_rim` | `S_s_dep` | `S_r_cri` | net demand | clamp | ground |
|---|---|---|---|---|---|---|---|---|
| mm/a | 538 | **8.1e+09** | **1.1e+10** | -7.5e+08 | **5.5e+10** | **7.35e+10** | **-7.35e+10** | **1166** |

**THE CLAMP REMOVES 99.999998 % OF THE DEMAND AND THE GROUND FLUX IS WHAT IS LEFT OVER.**
`P_snow` = 1046 mm/a is not the sum of any set of microphysical rates; it is `P_max_flux` and the
zero floor, applied level by level to a demand seven orders of magnitude larger. That is why
cutting the snow riming by 5 bought 62 %, the graupel riming by 5 bought 37 %, and a further 4.2
bought 11 %: **a saturated term does not care about its coefficient.** Every species pegging the
same 9.46e+04 mm/a per-cell ceiling — recorded above as "a cap binds somewhere in all three" — is
this, and the cap is `IceSchemeCommon::P_max_flux` = 3.0e-3 kg/(m2 s).

**THE CAUSE IS DIMENSIONAL: ThreeCat DIVIDES EVERY FLUX BY ITS OWN GROUND VALUE BEFORE USING IT.**

    Rain = P_rain[i] / max(P_rain[0], 1e-6)

and then hands that RATIO to relations that are dimensional in kg/(m2 s). The clearest case is
the mass content of falling rain, `r_q_r = A_r*B_rad^(-8/9)*Rain^(8/9)`: with the raw flux a
typical 1e-5 kg/(m2 s) gives **7.8e-5 kg/m3**, an ordinary 0.08 g/m3 of rain water; with the
ratio, which the 1e-6 floor lets reach `P_max_flux`/1e-6 = **3000**, the same expression returns
**~2.5e+3 kg/m3 — denser than liquid water.** Every collection term is a power of these
quantities, so all of them inflate together, and `S_r_cri ∝ ice/m_i * Rain^(13/9)` inflates most.

**TWOCAT USES THE RAW FLUX** (`TwoCatIceScheme.h:241`), and so does the COSMO documentation the
constants come from. The normalisation is ThreeCat's alone and has been there since the initial
commit. It is also the origin of the "ThreeCat NaN blow-up": `P_snow_0` in a denominator with
`S_s_rim ∝ Snow` is a geometric amplifier down the column, which is exactly what `P_norm_floor`
and `P_max_flux` were added to contain. **Both contain the symptom, and the containment is now
the model's precipitation.**

**AND THERE ARE NO AVAILABILITY LIMITERS AT ALL.** TwoCat carries five — a proportional
cloud-water limiter, a proportional ice limiter, and per-species clips on `S_ev`, `S_s_dep` and
`S_s_melt`. ThreeCat has none, and its own riming comments cite "the proportional cloud-water
limiter" as the mechanism they act through: **the comments were ported and the limiter they
describe was not.** Nothing stops a sink removing more water than the cell holds; the only bound
is the floor and the cap, applied one level later to the integrated flux. Measured on the
dimensionally-repaired branch, rain evaporation demands **5074 mm/a against 1664 mm/a of rain
sources** — three times the water present.

**THREE CHANGES, MEASURED SEPARATELY SO NONE IS CONFOUNDED** (same checkpoint, 6 iterations,
24 threads, mm/a):

| arm | P_rain | P_snow | P_graupel | total | vs NASA |
|---|---|---|---|---|---|
| shipped ThreeCat | 1118 | 1046 | 515 | **2679** | 2.74x |
| **+ in-cloud autoconversion (NOW THE DEFAULT)** | 661 | 1050 | 521 | **2081** | 2.13x |
| + `ATM_ICE_LIMITERS=1` | 52 | 22.6 | 1.5 | **76** | 0.08x |
| `ATM_ICE_RAW_FLUX=1` alone | 0.003 | 0.69 | 0.0 | **0.73** | 0.0007x |
| `ATM_ICE_RAW_FLUX` + in-cloud autoconversion | 93.0 | 0.69 | 0.0 | **93** | 0.10x |
| + `ATM_ICE_LIMITERS` (all three) | 8.7 | 0.69 | 0.0 | **9.4** | 0.010x |
| *TwoCat, same checkpoint, same 6 iterations* | *359* | *2.4* | — | *362* | *0.37x* |
| *NASA* | | | | *978* | |

**THE IN-CLOUD AUTOCONVERSION IS FLIPPED ON** — it is `3ea78e3`'s repair, which TwoCat, OneCat
and ZeroCat all received and ThreeCat did not (gap 2 above), and it is the same defect in the
same place: `c_c_au*(q_c - 0.0002)` tests a 0.2 g/kg threshold against the GRID MEAN, which under
the default `ATM_CLOUD_FRAC` peaks near 0.06, so ThreeCat's autoconversion was **off**. With the
fractional form `f*c_c_au*(q_c/f - q_c_crit)` and the matching accretion it takes total
precipitation 2679 -> 2081 and rain 1118 -> 661.

**THE OTHER TWO ARE DEFAULT OFF, AND THE REASON IS NOT DOUBT ABOUT THE PHYSICS.** A ratio fed to
a dimensional law is wrong and a sink that removes water that is not there is wrong. But each,
alone or together, takes precipitation to **10-90 mm/a** — the flip pattern this file has now
recorded six times, where removing one error of a cancelling pair makes the model worse. Read
them as measured, not as ready.

**WHAT THE REPAIRED BRANCH IS LIMITED BY, WHICH IS THE NEXT THING.** With the fluxes physical the
budget is legible for the first time: rain sources are `S_c_au` **1602** and `S_ac` 98, and the
sink is `S_ev` at **-5074** — evaporation removes essentially the whole flux at every level,
because a small flux evaporates faster than it falls (`S_ev ∝ Rain^(4/9)` against an availability
∝ `Rain`). TwoCat from the same state has the same shape and survives it: sources 2955, `S_ev`
demand 11562, **11.8 % surviving** against ThreeCat's 0.5 %.

**THE CLAMP IS NOW SPLIT INTO FLOOR / CAP / WINDOW, AND IT IS A DIFFERENT MECHANISM ON EACH
BRANCH AND FOR EACH SPECIES** (2026-09-01, same probe, mm/a):

| | net demand | floor (injects) | cap (truncates) | window (deletes) | ground |
|---|---|---|---|---|---|
| **shipped** snow | 2.174e+10 | +1.23e+07 | **-2.175e+10** | -1.26e+05 | 1166 |
| shipped graupel | 6.290e+09 | +0.0 | **-6.289e+09** | -7.51e+04 | 517 |
| shipped rain | -3.827e+09 | **+6.294e+09** | -2.466e+09 | **0.0** | 730 |
| **repaired** snow | 491.1 | +4.6 | **0.0** | **-499.0** | 0.7 |
| repaired rain | -5920.9 | **+5929.9** | **0.0** | **0.0** | 9.0 |

**ON THE SHIPPED BRANCH IT IS THE CAP, AND NOTHING ELSE.** 99.94 % of the snow removal and
essentially all of the graupel removal is `P_max_flux`; the window accounts for 0.0006 %. So
"`P_snow` = 1046 mm/a is `P_max_flux` applied level by level" is now measured rather than
inferred.

**AND THE SHIPPED RAIN IS NOT TRUNCATED, IT IS MANUFACTURED.** Its net demand is NEGATIVE —
-3.83e+09 mm/a, because `S_r_cri` removes rain at 6.15e+09 — and the floor **injects
+6.29e+09 mm/a of rain that no source produced**, six million times NASA's entire precipitation,
to keep the flux non-negative. The cap then removes 2.47e+09 of it. The rain reaching the ground
on the shipped branch is the residue of an injection, not of a generation.

**THE CAP IS INERT ONCE THE FLUXES ARE PHYSICAL** — 0.0 for every species on the repaired branch,
which is the independent confirmation that `ATM_ICE_RAW_FLUX` is the right diagnosis: no cell
reaches 3.0e-3 kg/(m2 s) any more.

**THE WINDOW IS THE WHOLE OF THE SNOW LOSS ON THE REPAIRED BRANCH, AND ZERO FOR RAIN — WHICH
HALF-REFUTES THE PARAGRAPH THIS ONE REPLACES.** Snow: the column generates 491 mm/a and the
window deletes **499**, leaving 0.7 at the ground. Snow crossing the melting level is discarded
rather than melted, in a scheme that HAS `S_s_melt`. That half stands and is now measured.
**The rain half was wrong**: the window removes **exactly 0.0**, because `S_c_au` and `S_ac` are
both gated on `t_u >= m.t_0`, so no rain is ever FORMED above the freezing level for the window
to delete. "Rain formed above the freezing level never reaches the ground" described a flux that
does not exist.

**WHAT LIMITED THE REPAIRED RAIN WAS THE FLOOR, BECAUSE THE LIMITER WAS CHARGING THE WRONG
QUANTITY.** The first version of the clip read `m.P_rain.x[i][j][k]`, which at that point in the
loop is the PREVIOUS pass's value at this level. The sink computed at level i is not charged
against that at all: it is charged one level lower, as
`P_x[i-1] = clamp(P_x[i] + rho[i]*S_x[i]*dz[i-1])`, so the flux it can take from is `P_x[i]` —
which the iteration is about to compute from `P_x[i+1]` and `S_x[i+1]`, both already known — over
a layer mass `rho[i]*dz[i-1]`, not `rho[i]*dz[i]`.

**AND THE FIRST REPAIR OF IT OVER-CLIPPED, BY A FACTOR OF TEN.** Charging a sink against the
ARRIVING flux alone forbids a level from evaporating the rain it makes in that same level, and
rain that condenses and evaporates inside one layer is a real process. The mass constraint is
only that the flux must not go NEGATIVE, and the level's own sources enter the same expression
the sinks are subtracted from, so the availability is **arriving + this level's sources**.
ThreeCat + `ATM_ICE_RAW_FLUX` + `ATM_ICE_LIMITERS`, same probe:

| clip charged against | `S_ev` | floor | cap | window | ground rain | total precip |
|---|---|---|---|---|---|---|
| stale same-level value (first try) | -7620 | **+5930** | 0.0 | 0.0 | 9.0 | 9.4 |
| arriving flux only (over-clips) | -1653 | +0.0 | 0.0 | 0.0 | 5.3 | 5.5 |
| **arriving + local sources** | **-1585** | **+0.0** | **0.0** | 0.0 | **72.0** | **72.9** |

**EVERY GRAM NOW LEAVES THROUGH A NAMED TERM.** Floor and cap are EXACTLY zero for all three
species; the only non-rate remover left in the scheme is the SNOW temperature window, -498.6 of
a generated 499.6.

### The accepted configuration's 992 mm/a stands on 8116 mm/a of manufactured water

**THE SAME DEFECT IS IN TWOCAT, THE DEFAULT SCHEME, AND IT IS MEASURED ON THE ACCEPTED
CONFIGURATION ITSELF** (`ATM_ICE_LIMIT_ARRIVING`, new, default 0 = shipped; `config_accept.xml`,
`nm` = 100 from scratch, moist gate 0, 24 threads — the run that produced the 992 mm/a every
cloud, humidity and radiation constant accepted on 2026-08-31 was fitted against):

| TwoCat, iteration 100 | sources | `S_ev` | `S_r_frz` | clamp | `P_rain(ground)` | **Precip mean** |
|---|---|---|---|---|---|---|
| **shipped** | 2846 | **7749 (272 %)** | 2231 | **+8116 INJECTED** | 366 | **992.2** |
| **arriving + local sources** | 3097 | **2347 (76 %)** | **388** | **-0.66 (0.0 %)** | 122 | **369.2** |
| *NASA* | | | | | | *978* |

**THE HEADLINE NUMBER DOES NOT SURVIVE MASS CONSERVATION, AND THE ARITHMETIC IS WORSE THAN A
PERCENTAGE.** `S_ev` demands **272 % of every rain source in the model**, so the shipped rain
budget is NEGATIVE — 2846 of sources against 9979 of sinks — and what reaches the ground is not a
reduced version of the sources, it is what `max(0, ...)` manufactured: **8116 mm/a, 8.2x the
reported total and eight times NASA's entire precipitation.** The identity closes on it:
2846 - 9979 + 8116 = **983**, against a reported `Precip mean` of 992.

*An earlier version of this paragraph called the total "63 % floor injection", which conflates
two numbers. 63 % is what is LOST when conservation is enforced (992 -> 369). The injection
itself is 8116 — larger than the reported total, not a fraction of it.*

With each sink able to take only water that is there, the same configuration rains **369 mm/a,
38 % of NASA**, and the clamp is **0.0 %**: the budget closes on rates alone for the first time
in this tree.

**THE OFF-BRANCH IS UNCHANGED, VERIFIED ON THE SAME RUN.** The control re-run with the new binary
reproduces the old one to every printed digit — sources 2846.40 against 2846.41, `Precip mean`
9.922e+02 against 9.922e+02, `P_rain mean` 9.704e+02 against 9.703e+02 — the documented
fixed-thread-count non-determinism and nothing else.

**SO THE PRECIPITATION AGREEMENT WAS A SIXTH CANCELLING PAIR, AND IT IS THE ONE THAT MATTERED
MOST.** An over-large evaporation sink and a floor that manufactured water to feed it. Both of
this tree's NASA matches — the shipped 1046 and the accepted 992 — are now explained, and neither
was a result. **Do not quote either as validation.** What can be said after today is narrower and
firmer: the microphysics can be made mass-conserving, and when it is, this model rains about a
third of what Earth does.

### The evaporation was un-weighted, and weighting it recovers the whole number

**`S_ev` IS THE SINK THAT FORCES THE FLOOR, AND IT WAS SPREAD OVER THE WHOLE GRID BOX**
(`ATM_RAIN_AREA=<fraction>`, new, default 0 = shipped, both schemes). Rain falls in shafts, not
evenly; evaporation is local, so the grid-mean tendency of a shaft covering a fraction `f` is
`f*E(R/f)`, and because `E ∝ R^(4/9)` is SUB-linear the transform reduces it by `f^(5/9)`.
**Fifth occurrence of the grid-mean defect** — the same shape as the autoconversion (`3ea78e3`),
the saturation adjustment (`6d163a0`), the convective trigger and the updraft condensate. The
`(q_sat - c)` driver is left grid-mean, because the in-shaft humidity is not a quantity this
model carries, so the correction is a LOWER bound.

`config_accept.xml`, `nm` = 100 from scratch, 24 threads, mm/a:

| arm | sources | `S_ev` | clamp | **floor injected** | **Precip mean** |
|---|---|---|---|---|---|
| **shipped** | 2846 | **7746 (272 %)** | -8114 | **8128 (819 % of the total)** | **992.8** |
| shipped + area 0.10 | 2825 | 2315 (82 %) | -2445 | 2458 (279 %) | 880.5 |
| mass-conserving, no area | 3097 | 2347 (76 %) | -0.66 | ~0 | **369.2** |
| + area 0.30 | 2770 | 1493 (54 %) | -0.03 | 0.44 | 892.9 |
| **+ area 0.10** | 2718 | **1362 (50 %)** | **-0.01** | **0.20** | **975.8** |
| + area 0.05 | 2749 | 1494 (54 %) | 0.00 | 0.10 | 877.6 |
| *NASA* | | | | | *978* |

**A MASS-CONSERVING CONFIGURATION LANDS ON NASA.** `ATM_ICE_LIMIT_ARRIVING=1 ATM_RAIN_AREA=0.10`
gives **975.8 mm/a against 978, with the floor injecting 0.20 mm/a instead of 8128** — the first
time this model's headline output has been a difference of rates rather than a residue of
manufactured water. The identity closes exactly: 2718 - 1749 + 0.01 = **968.46**, the printed
`P_rain(ground)` to the digit.

**THE TWO KNOBS ARE A PAIR AND EITHER ALONE IS WORSE.** Conservation alone gives 369; area
weighting alone leaves 2458 mm/a of injection standing and still moves the total to 880. Sixth
time in this file that fixing one of a cancelling pair makes the model worse.

**AND THE AGREEMENT SITS ON A PEAK, NOT A PLATEAU — READ IT AS A FITTED CONSTANT.** 877.6 at
0.05, 975.8 at 0.10, 892.9 at 0.30: the response is NON-MONOTONIC and one step either way costs
10 %. Two effects oppose each other — a smaller `f` cuts `S_ev` by `f^(5/9)` but raises the flux
density `R/f` that `S_ev` is computed FROM, and the larger surviving flux then feeds back through
accretion (55 -> 84 -> 125 as `f` falls) and snowmelt (89 -> 18 -> 7.5). **0.10 was chosen
because it lands on 978**, and a single constant cannot represent both stratiform rain, which
covers most of a 1-degree box, and convective rain, which covers a few per cent. This is a
scaffold in exactly the sense `ATM_RH_MIN` is one.

**FLIPPED ON BY DEFAULT 2026-09-01, AT THE USER'S INSTRUCTION, AS A PAIR.**
`ATM_ICE_LIMIT_ARRIVING` 0 -> 1 and `ATM_RAIN_AREA` 0 -> 0.10. Both directions verified on
`config_accept.xml`, 24 threads:

| | `S_ev` | clamp | floor injected | `Precip mean` |
|---|---|---|---|---|
| **new default**, clean environment | 1362 (50.1 %) | -0.01 | **0.199 (0.0 %)** | **975.8** |
| `ATM_ICE_LIMIT_ARRIVING=0 ATM_RAIN_AREA=0` | 7748 (272 %) | -8116 | 8129 (819.2 %) | **992.4** |

The new default reproduces the measured arm to every printed digit (975.8, `P_rain(ground)`
968.46, `S_ev` 1362.00) and the restore reproduces the shipped branch to 0.04 % — the documented
fixed-thread-count non-determinism. **Setting those two variables back restores the old branch**,
and the `[RUN CONFIG]` banner prints both with `*` when they are compiled-in.

### Scored against the NASA FIELD at last, and the mean was hiding a tropical spike

**THE MODEL HAS BEEN COMPARED WITH `precipitation_NASA` BY ITS GLOBAL MEAN ALONE SINCE THE
BEGINNING.** The field is read into memory, written to VTK and printed as one scalar beside the
model's; the two have never been differenced. `printDataAtm` now scores them — cos-lat weighted
on the same weights `GetMean_3D` uses, so the model mean it reports IS the `Precip mean` above
it. Unconditional, one 2-D pass, `nm` = 100:

| | model mean | bias | **pattern r** | centred RMS | sigma model/NASA |
|---|---|---|---|---|---|
| **new default** | 975.8 | **-0.3 %** | **+0.456** | **1432** | **2.33** |
| `LIMIT_ARRIVING=0 RAIN_AREA=0` | 993.6 | +1.6 % | **+0.214** | 2016 | 2.97 |
| *NASA* | *978.3* | | | | |

**THE FLIP IS VINDICATED BY SOMETHING IT WAS NOT TUNED AGAINST.** The pattern correlation more
than DOUBLES, 0.214 -> 0.456, the centred RMS falls 29 % and the variance ratio moves toward 1 —
and none of that was visible in a mean that was already "right" on both branches. The land/ocean
partition is the clearest single line:

| mm/a | model, new default | model, old branch | NASA |
|---|---|---|---|
| **land** | **766.3** | **2217.0** | **782.3** |
| **ocean** | **1058.7** | **509.4** | **1055.8** |

**The old branch was 2.8x too wet over land and 2x too dry over ocean, and its global mean was
right to 1.6 %.** The injection was concentrated where the sub-cloud air is driest — over land —
so removing it repaired the partition almost exactly. That is the strongest evidence the flip is
a physics repair and not a re-tune: nothing in it knew about the land mask.

**AND THE REMAINING ERROR IS ALL ONE SHAPE: A TROPICAL SPIKE WITH NOTHING POLEWARD OF IT.**

| mm/a by \|latitude\| | 0-15 | 15-35 | 35-65 | 65-90 |
|---|---|---|---|---|
| **model** (new default) | **3352.5** | **261.9** | **157.1** | **7.6** |
| model (old branch) | 2611.5 | 954.6 | 111.1 | 7.2 |
| *NASA* | *1487.0* | *761.4* | *981.1* | *364.2* |

2.3x too much rain in the deep tropics, 3x too little in the subtropics, **6x too little in the
storm-track band and 48x too little at the poles.** `sigma` = 2.33 is that spike measured a
second way. **This is a CIRCULATION result, not a microphysics one**: mid-latitude precipitation
is made by baroclinic eddies, and this tree's recorded jet spin-down, its `Psi(ground)`
non-closure at ~40 % of the Hadley cell, and its dead convective scheme are exactly the items
that would starve 35-65 degrees. The microphysics work of 2026-09-01 has taken the water budget
as far as it can go on its own.

**SO THE GLOBAL MEAN WAS NEVER THE PROBLEM AND IS NOT THE TARGET.** Two configurations match it
to 2 % with pattern correlations of 0.21 and 0.46 and land/ocean splits that differ by a factor
of three. Quote `r`, `sigma` and the four bands, or quote nothing.

**WHAT IT IS AND IS NOT.** It IS: the microphysics conserving mass, one cancelling pair removed,
and the number no longer standing on an injection. It is NOT a validation — `nm` = 100 is 20
seconds of physical time, the 1000-iteration run had precipitation still climbing, and the
PATTERN has never been compared with the NASA field at all, only its global mean.

**THE FLOOR INJECTION IS NOW PRINTED IN EVERY RUN**, unconditionally, directly under
`Precip mean`, as an absolute figure and as a percentage of the reported total. On the shipped
branch it reads 819 %. It is a function-local static in `IceSchemeCommon` rather than a model
member, deliberately: adding a member moves `sizeof(cAtmosphereModel)` and that is the
stack-canary hazard.

**AND `ATM_SR_DIAG`'s GROUND BUCKET WAS UNDER-COUNTING, WHICH IS HOW THE IDENTITY EXPOSED IT.**
The budget closed on 983 while the print's own `P_rain(ground)` read 366. The ground value was
assigned AFTER the convergence test, and that test `break`s — so every column that settled on
that pass skipped its own assignment and contributed zero. Fixed: it is written before the test.
**The sources, sinks and clamp were never affected**, so every number quoted from this instrument
stands; only its ground column was low. Fourth instrument-shaped defect found by making a budget
close.

**AND THE FIX DOES NOT CLOSE TWOCAT COMPLETELY WHERE IT CLOSES THREECAT.** `S_r_frz` falls
2231 -> 388 because it is now clipped with the rest, but TwoCat's shipped code has no clip on it
at all; the residual -0.66 is 0.02 % and the remaining gap to NASA is a physics gap, not a clamp.

**READ EVERY NUMBER IN THIS SUBSECTION AS A SIX-ITERATION RESPONSE TO A FIELD THAT WAS SPUN UP BY
THE UNREPAIRED SCHEME**, not as a climate. The checkpoint's cloud, ice and vapour were made by
1e+10 mm/a of demand being clamped; six iterations is ~1.2 seconds of physical time and cannot
re-equilibrate them. The arms are comparable to each other because they start from one state, and
none of them is a prediction.

**SO THE DEFAULT WENT BACK TO TWOCAT ON 2026-09-01**, at the user's instruction, on this
evidence: the calibration record through 2026-08-31 was all measured against TwoCat, and TwoCat
is the only scheme here whose precipitation is made of rates rather than of caps. `= 3` restores
ThreeCat.

**AND TWOCAT HAS THE SAME DISEASE, FAR MILDER, WHICH THE REVERT DOES NOT CURE.** From the same
checkpoint (`ATM_SR_DIAG`, 6 iterations): sources 2955, `S_ev` demand **11562 -- 391 % of the
sources** -- clamp -11816, ground 350, **11.8 % surviving** against ThreeCat's 0.5 %. So
TwoCat's answer is also partly set by the `max(0, ...)` floor, by a factor of ~4 rather than
~1e+7, and its `S_ev` carries no rain-area fraction either. Reverting buys a scheme whose numbers
are MOSTLY rates, not one that is clean.

## The jet spin-down is the radial Shapiro filter, measured two ways and 100 % attributed

**THE JET IS NOT SPINNING DOWN. IT IS BEING FILTERED AWAY** (2026-09-01). Over the 1000-iteration
run the zonal-wind maximum falls **27.58 -> 18.55 m/s** and its core **sinks 9007 -> 5512 m**, in
a quasi-linear decay that shows no sign of equilibrating.

**THE ATTRIBUTION WAS ALREADY WRITTEN AND NOBODY HAD OPENED IT.** That run left 45
`w_momentum_budget_*.csv` files in `output_spin1000/`. At the jet core (31S, 4988 m,
`wbar` = 16.48 m/s), per iteration:

| stage | dw, m/s |
|---|---|
| **`dw_radial`** — the radial (vertical) Shapiro filter | **-3.707e-03** |
| `dw_dyn` — ALL the physics: PGF + Coriolis + advection + diffusion + drag | **-8.47e-06** |
| `dw_polar`, `dw_orog` | 0.000 — the jet sits at 31 deg, outside both |
| sum of the four captured stages | -3.715e-03 |
| **actual `wbar` change, 980 -> 1000** | **-3.750e-03** |

The budget closes to 99.1 %, and a NUMERICAL SMOOTHER is **99.8 % of what it captures, 437x the
entire resolved dynamics.** *Seventh occurrence of "the answer was in output nobody opened".*

**AND THE A/B CONFIRMS IT EXACTLY, ON A KNOB THAT WAS ADDED FOR THIS ON 2026-07-21 AND NEVER
SWEPT.** `ATM_RADIAL_SHAPIRO_STRENGTH` (default 1.0) appears in no results table in this file.
Four arms, 600 -> 700 from the same checkpoint, 24 threads:

| strength | jet at 620 | jet at 700 | decay over 80 iters | per 20 iters |
|---|---|---|---|---|
| **1.0 (shipped)** | 18.4293 | 17.9733 | **0.4560** | 0.1140 |
| 0.5 | 18.4842 | 18.2586 | 0.2256 | 0.0564 |
| 0.25 | 18.5120 | 18.3941 | 0.1179 | 0.0295 |
| **0.0 (off)** | 18.5396 | **18.5388** | **0.0008** | 0.0002 |

**EXACTLY LINEAR IN THE STRENGTH** — 0.456 x 0.5 = 0.228, x 0.25 = 0.114, x 0 = 0 — and with the
filter off the jet is FROZEN to four decimal places, 570x less decay. **The spin-down is the
filter and nothing else.**

**THE MECHANISM IS THE STRETCHED GRID, AND IT EXPLAINS THE SINKING CORE.** `dw_radial` down the
31S column: **+9.1e-04 at 1211 m, +2.0e-03 at 2163 m** — the filter ADDS momentum on the lower
flank — then **-3.7e-03 at the 4988 m core and -5.6e-03 at 6719 m**, removing most ABOVE the
maximum. A smoother eating a vertical maximum takes from the peak and gives to the flanks; the
removal is top-heavy because the filter is uniform in INDEX space while the layer spacing runs
from 38.9 m at the bottom to over a kilometre aloft. **This tree's stretch spread is 23.21x, the
worst in the family**, so the same filter is far harsher in physical terms at the top of the jet
than at its base — and the maximum migrates downward, 9007 -> 5512 m.

**AND THERE IS NO GENERATION TERM, WHICH THE FILTER HAS BEEN HIDING.** At strength 0 the jet does
not recover, it FREEZES: nothing damps it and nothing drives it. **Removing the filter stops the
decay; it does not give the model a jet that is maintained by a temperature gradient**, and the
thermal-wind balance that maintains a real jet is not acting here.
**THE `pgf` = 1.5e-11 THAT WAS QUOTED HERE AS THE EVIDENCE IS AN IDENTITY, NOT A MEASUREMENT**
(2026-09-02) — it is the ZONAL pressure gradient in a ZONAL MEAN, which is zero over a periodic
domain by construction. The CONCLUSION survives and is now measured in the right equation: see
*There is no thermal wind* below, where the MERIDIONAL `pgf`/`coriolis` is 2.5e-05 and the cause
is that no term in the horizontal momentum equations carries the temperature at all.

**THE FROM-SCRATCH TEST HAS NOW BEEN RUN, THE GUARD HOLDS AT QUARTER STRENGTH, AND THE JET IS
23 % STRONGER — AND IT BUYS NOTHING IN MID-LATITUDE RAIN** (2026-09-01, default TwoCat
configuration, `config_accept.xml` at `nm` = 600, moist physics from iteration 0, 24 threads,
62 min per arm). Zonal-mean jet above 3 km, from the run's own budget files:

| iteration | 20 | 160 | 320 | **600** | core height at 600 |
|---|---|---|---|---|---|
| **strength 1.0** | 26.488 | 22.573 | 20.422 | **18.303** | **5513 m** |
| **strength 0.25** | 27.190 | 25.590 | 24.277 | **22.599** | **6719 m** |

**BOTH ARMS RAN 600 ITERATIONS FROM SCRATCH, EXIT 0, ZERO NaN** — past 155, 357 AND 483, the
three documented failure points. `max|u|` ends at 0.0412 against the control's 0.0237: larger, as
expected from a weaker filter, and three orders below anything resembling the +-100 m/s coastal
blow-ups the guard exists for. **At quarter strength the guard still holds.**

**BUT THE ATTRIBUTION IS WEAKER HERE THAN THE RESTART TEST SAID, AND THE TWO DISAGREE.** Over the
settled late stretch (iterations 400 -> 600) the decay is **0.674 per 100 iterations at strength
1.0 and 0.517 at 0.25**. Fitting `decay = A*s + B`, the filter-proportional part is A = 0.210 and
the residual **B = 0.464 — so the filter is only 31 % of the late decay here**, against ~100 % in
the 600 -> 700 restart test where strength 0 froze the jet outright. The two measurements are on
different configurations (ThreeCat restart vs TwoCat from scratch) and different fields, and the
from-scratch core is still migrating downward in BOTH arms at iteration 600, so part of `B` is
plausibly the analytic initial condition still adjusting rather than a steady drain. **The
discrepancy is real and unexplained; the clean resolution is a strength-0 arm on THIS
configuration, and it has not been run.**

**AND THE STORM TRACK DOES NOT COME BACK. AT ALL.**

| iteration 600 | Precip | pattern r | 0-15 | 15-35 | **35-65** | 65-90 |
|---|---|---|---|---|---|---|
| strength 1.0 | 998.4 | +0.457 | 3432.2 | 271.0 | **156.6** | 7.5 |
| strength 0.25 | 1001.4 | +0.461 | 3443.0 | 272.0 | **156.6** | 7.5 |
| *NASA* | *978.3* | | *1487.0* | *761.4* | ***981.1*** | *364.2* |

**156.6 in both arms, identical to four figures**, with `r` moving 0.457 -> 0.461 and the
land/ocean split unchanged. **A 23 % stronger, 1.2 km higher jet buys nothing.** That is the
prediction made before the run and it is confirmed: the jet is unmaintained whether or not it is
being eroded, and **the eddies are STARVED, not damped.** (The `pgf` = 1.5e-11 first cited here is
the wrong instrument — see below — but the eddy statement stands on `wbud_advh` = +4.0e-08 against
`wbud_cor` = -3.2e-06, 80x smaller.)
Easing the filter preserves the shear; it does not create the baroclinic growth that turns shear
into rain. **Do not expect a filter setting to fix the 35-65 degree band.**
**AND THE 156.6 IN BOTH COLUMNS IS NOT EVIDENCE ABOUT THE EDDIES AT ALL — corrected 2026-09-03.**
`ATM_TW_BALANCE` has since produced a jet **71 % stronger and relocated INTO this band**, and the
35-65 deg precipitation moved **0.013 %**. 600 iterations is 120 s, over which a 0.33 m/s meridional
flow carries a parcel **36 metres**, so no dynamical arm of this length can move a precipitation
band whatever it does to the circulation. The eddy statement stands where it was measured, in the
momentum budget; the precipitation columns of this table measure the advective displacement. See
*The balanced initial state DOES give this model a mid-latitude jet* below.

*Useful by-product*: the control is `Precip` 998.4 at `r` = 0.457 after 600 iterations against
975.8 at `r` = 0.456 at `nm` = 100 — so under the 2026-09-01 defaults the precipitation field is
STABLE between iterations 100 and 600, where the pre-flip 1000-iteration run climbed 992 -> 1215.
And both arms wrote `atm_restart_0Ma_600.bin`, the first spun-up checkpoints on the CURRENT
default configuration; every probe before this restarted from a ThreeCat-conditioned state.

**AND THE PRECIPITATION NUMBERS IN THOSE FOUR ARMS ARE WORTHLESS — WRONG CONFIG, MY ERROR.** They
were run from `config_ggrim.xml`, which is a ThreeCat probe with the raw-flux and limiter knobs
OFF, so all four report ~13300 mm/a at `r` = -0.02: ThreeCat's clamp residual, identical across
arms because the clamp does not care about the jet. **The storm-track question — does easing the
filter recover the 35-65 deg band, 157 against NASA's 981 — is UNANSWERED**, needs the default
TwoCat configuration, and needs far longer than 100 iterations because baroclinic eddies have to
grow before they rain.

## There is no thermal wind, because nothing in the horizontal momentum equations carries the temperature

**THE QUESTION THIS FILE ASKED — "why is `pgf` = 1.5e-11 at the jet core?" — HAS AN ARITHMETIC
ANSWER, AND THE PHYSICS QUESTION IS NEXT DOOR** (2026-09-02).

`wbud_pgf` is `-dpdphi_invrs`, the ZONAL pressure gradient, and `write_w_momentum_budget` reports
it as a ZONAL MEAN. The zonal mean of a zonal derivative over a periodic domain is zero to
round-off, so 1.5e-11 is an identity, not a measurement: median `|wbud_pgf|/|wbud_cor|` over
20-70 deg above 3 km is **4.2e-06**. **The zonal-mean zonal momentum equation has no
pressure-gradient term by construction, and nothing about the jet's forcing can be read off that
print.** What CAN be read off it is that the jet's two real sources are Coriolis on the mean
meridional flow (-3.16e-06) and eddy momentum flux convergence (`adv_horiz`, **+3.98e-08, 80x
smaller**) — which is the "eddies are starved" statement, measured.

**THE BALANCE THAT MAINTAINS A JET IS THERMAL WIND, AND IT LIVES IN THE MERIDIONAL EQUATION.
THERE IS NONE.** From the same run's `v_momentum_budget_600.csv` (`output_rsfull100`, the
from-scratch default-configuration checkpoint), 20-70 deg, above 3 km, n = 1632:

| v-equation term | median \|term\| | at the jet core (31S, 5513 m) |
|---|---|---|
| `coriolis` | 1.390e-04 | -2.7392e-04 |
| `diffusion` | 6.500e-08 | +4.78e-08 |
| `adv_vert` + `adv_horiz` | ~8e-09 | -2.6e-08 |
| **`pgf`** | **2.769e-09** | **+2.7413e-08** |
| **\|pgf\|/\|coriolis\|** | **2.49e-05** (p95 3.2e-03) | **1.00e-04** |

**The meridional momentum equation is `dv/dt = -f w` and nothing else.** Coriolis is 50 000x the
pressure gradient. The model is not geostrophic anywhere, at any latitude, at any height.

**AND THE MODEL'S OWN TEMPERATURE FIELD ASKS FOR SIX TIMES THE JET IT HAS.** At 31S on the 87E
section, `dT/dy` = **0.683 K per degree** over 1-10 km, which by thermal wind
(`du/dz = -(g/fT) dT/dy`) calls for **+28.1 m/s** of westerly shear over that depth. The model's
own shear there is **+4.56 m/s — 16 %**. The shear is present in the temperature and absent from
the wind.

### The cause is one sentence, and it is structural rather than a mis-scaling

**NOTHING IN `rhs_u`/`rhs_v`/`rhs_w` CARRIES THE TEMPERATURE FIELD.**

1. **The only pressure in the momentum equations is `p_dyn`** (`RHS_Atm_Turb.cpp:95-340` — every
   one of the twelve `dpdr`/`dpdthe`/`dpdphi` branches reads `p_dyn` and nothing else).
   **`p_stat` — which `ThermoAtm` builds barometrically and which holds the entire hydrostatic,
   thermal structure — appears in no momentum equation.**
2. **`p_dyn` is a divergence-removal pressure, not a mass-field pressure.** Its Poisson source is
   `div(aux)`, the divergence of the momentum TENDENCY (`PressureSolverAtm.h:427`). A projection
   pressure removes the divergent part of the forces already present; it cannot manufacture the
   geostrophic pressure a mass field implies, because no mass field is given to it.
3. **The Boussinesq buoyancy is the one surrogate route from `t` into momentum, it is RADIAL
   only, and it is inert.** Measured in this very run (`output_pmf_ctl`, iteration 700):
   `ubud_buoy` = **2.46e-05** against `ubud_pgf` = **8.23e-02** — **0.03 %, a factor of 3346.**
   `ATM_BUOY_CONSISTENT` puts the correctly non-dimensionalised value 5e5 higher (the 400 000x
   recorded in the knob section above).

So thermal-wind balance is impossible here **by construction**, and the jet is what
`VelocityInitializer`'s analytic profile left behind, being eroded by the radial Shapiro filter.
**That closes yesterday's result from the other side**: easing the filter bought a 23 % stronger,
1.2 km higher jet and moved the 35-65 deg band by 0.0 mm/a, because a jet with no thermal-wind
generation is a decaying initial condition whichever rate it decays at.

### The obvious suspect was swept, and it is a null — `ATM_POISSON_METRIC_FIX`

**THE POISSON OPERATOR'S HORIZONTAL METRIC REALLY IS INCONSISTENT WITH ITS OWN div AND grad.**
The source (`:427`) and the gradient correction (`:1019` in the time loop, `:1092` in
`project_initial_velocity`) each carry ONE power of `inv_rm`, so the composition `div.grad`
carries two; the Laplacian (`:240-241`) carries one. The RADIAL term already does it right —
`exp_2_rm` in the Laplacian against `exp_rm` once each in div and grad. **So the horizontal
pressure response is under-resolved by `rmet` = 397.5 by construction**, and the knob has existed
since 2026-07-21 and appeared in no results table. Swept at last, 600 -> 700 from
`output_rsfull100/atm_restart_0Ma_600.bin`, `config_accept.xml` defaults, 24 threads, exit 0:

| iteration 700 | control | `ATM_POISSON_METRIC_FIX=1` |
|---|---|---|
| v-eq `pgf` at the jet core | +2.850e-08 | **+3.339e-08 (+17 %)** |
| `pgf`/`coriolis` | 1.075e-04 | 1.259e-04 |
| **jet, m/s** | **17.720** | **17.720** |
| `div(u)` rms | 1.888e-03 | 1.888e-03 |
| `Psi(ground)` | — | identical to six figures |
| **Precip / r / all four bands / land / ocean** | 1014.9 / +0.458 | **identical to every printed digit** |
| rms `p_dyn` | 6.412e-04 | 6.560e-04 |
| `p_dyn` checkerboard index | 0.0164 | **0.0429 (2.6x)** |

**IT BUYS 17 % WHERE THE OPERATOR ALGEBRA PREDICTS 398x, AND THE REASON RETIRES THE KNOB.** The
Poisson operator is RADIALLY DOMINATED. At the jet level `c_r = 2*exp_2_rm/dr^2` = **431** against
`c_the = 2*inv_rm/dthe^2` = **16.5** — a ratio of 26 shipped, and **10 400 with the fix**. The
pressure amplitude is set by `c_r*k_r^2`, so the horizontal Laplacian coefficient is not the lever
on the horizontal pressure gradient at all; making it consistent moves the pressure by 2 % and the
velocity by nothing, while multiplying the checkerboard by 2.6. **Default stays 0.** (This is also
the resolution of `c_phi/c_r` = 0.0322 recorded above: the horizontal directions are weakly
constrained, and making them *more* weakly constrained is the consistent choice and changes
nothing.)

### And `p_dyn` is not in the units this tree has always said it is

**THERE IS NO EULER NUMBER IN THE PRESSURE-GRADIENT TERM.** `rhs_v` gets `-inv_rm * dp/dthe` with
no `p_0/(rho u_0^2)` coefficient, so for the momentum equation to be dimensionally consistent
`p_dyn` must be non-dimensionalised by **`rho*u_0^2` ~ 43.7 Pa**, not by `p_0` = 1013.25 hPa.
`PressureSolverAtm.h:98` states the wrong one, and **every VTK writes `p_dyn * p_0`, i.e. 2320x
too large.** Verified three independent ways:

| check | predicted | measured |
|---|---|---|
| the code's `pgf` rebuilt from the written `p_dyn` | 3.558e-05 | 3.427e-05 (zonal mean vs 87E) |
| the code's Coriolis rebuilt from `force_nd` | -0.3446 | -0.3424 |
| **`pgf`/`cor` from the field under `rho u_0^2`** | **1.03e-04** | **1.00e-04 (the budget's own)** |

Under the `p_0` normalisation the third line would read **0.240**, and the budget would have to be
wrong. It is not; the label is.

**CONSEQUENCE: THE TWO PRESSURE CLAMPS ARE THREE ORDERS SMALLER THAN THEIR COMMENTS SAY.**
`p_dyn_cap` = 2.0 is **87 Pa**, not 2000 hPa; `p_dyn_ceiling` = 3.0 is **131 Pa**. A
geostrophically balanced mid-latitude pressure field in these units is `p_dyn` ~ **70** — **23x
above the ceiling.** **Neither binds today** — max \|`p_dyn`\| is 0.017 at iteration 600, 176x
below the ceiling — so they are a LATENT barrier and not the present cause. But any repair that
gives this model a real thermal pressure field will hit them before it works, and re-sizing them
is part of that repair rather than a follow-up. *Sixth instrument-shaped defect in this tree,
after `Psi`'s constant density, `Q_Sensible`, `brunt_N2`'s terrain extrema, the "OLR" that was the
lid temperature, and `ATM_SR_DIAG`'s ground bucket.*

### The hydrostatic term is now written: `ATM_HYDRO_PGF`, default 0.0 and byte-identical off

**THE REPAIR NAMED ABOVE IS IN** (2026-09-02, `RHS_Atm_Turb.cpp`). A strength, not a flag:
`ATM_HYDRO_PGF=<s>` adds `-s*(100/(rho*u_0^2))*inv_rm*d(p_hyd)/dthe` to `rhs_v` and the matching
`inv_rmsinthe` term to `rhs_w`. **Off-branch verified BYTE-IDENTICAL over all 6 physics files at
1 thread**; the only file that differs is `RUN_CONFIG.txt`, by the new `[RUN CONFIG] dynamics
knobs` line (which also prints `POISSON_METRIC_FIX`, `RADIAL_SHAPIRO_STRENGTH`, `V_MASSBAL` and
`BUOY_CONSISTENT`, none of which any run was recording).

**IT CARRIES AN EULER NUMBER AND `p_dyn` DOES NOT, AND THAT IS THE WHOLE POINT OF THE SUBSECTION
ABOVE.** `p_dyn` is normalised by `rho*u_0^2` so its three `dpd*` lines need no coefficient;
`p_stat` is a real pressure in hPa, so this term needs `100/(rho u_0^2)` and gets it.

**HORIZONTAL ONLY, AND THE RADIAL OMISSION IS THE MEASURED PART, NOT AN OVERSIGHT.** In the
radial direction the hydrostatic gradient is gravity to within a per cent, and `rhs_u` does not
contain `-g` — the model is Boussinesq, so `rhs_u` carries only the RESIDUAL, as
`buoyancy*g*dt/u_0*(t - t_ref)`. Measured at 31S on the 87E section:

| z | `-(1/rho)dp/dz` | g | residual |
|---|---|---|---|
| 938 m | +9.7920 | 9.81 | **-0.018** |
| 5513 m | +9.7230 | 9.81 | **-0.087** |
| 9923 m | +9.6377 | 9.81 | **-0.172** |

In the model's own non-dimensional radial tendency (x `L/u_0^2` = 250.4) the full gradient is
**+2434.3** against `ubud_pgf` = 0.0823, `ubud_diff` = 0.00235 and `ubud_buoy` = 0.0000246 —
**29 586x the largest term in `rhs_u`**, which unbalanced would accelerate every cell upward by
**1.95 m/s per iteration** and peg the +-100 m/s clamp globally in ~50 iterations. The
normalisation does not save it either: dividing by `p_stat(1)` removes the horizontal offset, not
the vertical structure. **Adding it radially means adding `-g` with it, and their difference IS
the buoyancy term already there.** So the radial direction was the one place this model was
already doing the right thing. Do not "complete" the port with a third line.

**TWO FORMS, AND THE BAROCLINIC ONE IS THE DEFAULT WHEN THE KNOB IS ON.** `ATM_HYDRO_PGF_RAW=1`
differences `p_stat` as `ThermoAtm` builds it. That is the literal reading, and it is wrong at the
surface for a reason worth recording on its own: **`densities()` sets the sea-level pressure to
`p_sl = 1e-2*r_air*R_Air*T_surface`** — a CONSTANT-DENSITY sea-level pressure, correlating
**+1.0000** with surface temperature by construction and spanning **749..1046 hPa** against
Earth's ~980..1040. The default instead normalises every column at a common HEIGHT,
`p_hyd = p_0*p_stat(i)/p_stat(1)`, which leaves the BAROCLINIC part — the thermal wind. Level 1
rather than level 0 because `densities()` overwrites `p_stat.x[0]` with the value at
`i_topography`, and normalising a column at its own ground is the sigma-coordinate error (Tibet's
ground is at ~4 km: `p_hyd` at 5 km would be 0.88`p_0` against 0.51`p_0` over the neighbouring
ocean). Level 1 is never overwritten and holds the true barometric value in rock as in air.
**`p_hyd` is then an analytic function of (height, surface temperature) alone, so it is
horizontally smooth across every coastline and plateau edge and cannot produce the
sigma-coordinate pressure-gradient error at all.** The price is that it knows nothing about
topography: there is no mountain in this pressure field, only a temperature pattern.

**WHAT IT GIVES, READ OFF THE MODEL'S OWN FIELD BEFORE IT WAS RUN** (87E, 31S, dp/dy northward in
hPa/deg, and the zonal wind that gradient is in geostrophic balance with):

| z | 39 m | 938 | 2163 | 5513 | 9923 |
|---|---|---|---|---|---|
| raw `dp/dy` | 2.143 | 2.134 | 2.100 | 1.906 | 1.533 |
| **baroclinic** | **-0.000** | **0.208** | **0.438** | **0.816** | **0.937** |
| -> `u_g` | -0.00 | +2.27 | +5.41 | **+14.32** | **+26.88** |
| the model's own `w` | +0.62 | +2.51 | +7.29 | **+17.75** | **+6.60** |

**THE FLOW IS ALREADY CLOSE TO THERMAL-WIND BALANCE THROUGH THE LOWER AND MIDDLE TROPOSPHERE —
to ~20 % — AND NOTHING WAS HOLDING IT THERE.** That is the finding restated from the other side.
The exception is the top row: at 9.9 km the temperature calls for **26.9 m/s** and the model has
**6.6**, so this term builds the upper-tropospheric jet the 35-65 deg storm track needs.

**THREE ARMS, 600 -> 700 from `output_rsfull100/atm_restart_0Ma_600.bin`, `config_accept.xml`
defaults, 24 threads, all exit 0 with ZERO NaN** (`max|w|` 18.0-18.1, `max|u|` 0.022, i.e. the
CFL guard is nowhere near touched at FULL physical strength):

| iteration 700 | control | **`ATM_HYDRO_PGF=1.0`** | `+ _RAW=1` |
|---|---|---|---|
| v-eq `pgf` at the jet core | +2.85e-08 | **+1.989e-04** | +4.55e-04 |
| `pgf`/\|`coriolis`\| there | 0.0001 | **0.7502** | 1.7144 |
| **ageostrophic residual there** | **0.9999** | **0.2498** | 0.7144 |
| band median `pgf`/\|`cor`\| | 0.0000 | 2.209 | 4.378 |
| band median residual | 1.000 | 1.209 | 3.378 |
| **band p05 residual** | **0.997** | **0.079** | 0.615 |
| median \|`pgf`\| below 1 km | 3.04e-09 | **1.48e-05** | **3.43e-04** |
| jet, m/s | 17.7201 | **17.7201** | 17.7222 |
| Precip / r | 1014.9 / +0.458 | 1014.8 / +0.458 | 1014.9 / +0.458 |

**THE MODEL HAS GEOSTROPHICALLY BALANCED CELLS FOR THE FIRST TIME.** The p05 of the ageostrophic
residual goes **0.997 -> 0.079**: before, the MOST balanced cell in the extratropical free
troposphere was 99.7 % out of balance; now the balance closes to 8 %. At the jet core Coriolis is
75 % opposed where it was 0.01 % opposed.

**AND THE RAW BRANCH IS WORSE ON EVERY MEASURE, WHICH IS THE PREDICTION CONFIRMED.** Residual at
the core 0.714 against 0.250, band p05 0.615 against 0.079, and **23x the baroclinic term's
near-surface pressure gradient** — that last row is `p_sl ∝ T_surface` measured directly, a
barotropic surface gradient that should not exist. Use the default.

**WHAT IT DOES NOT DO, AND THIS IS THE LIMIT THAT MATTERS: THE VELOCITY DOES NOT MOVE.** The jet
is 17.7201 in both the control and the full-strength arm — identical to six figures — and every
precipitation number is identical to the printed digit. That is not a failure of the term, it is
the timescale recorded above: geostrophic adjustment takes `1/f`, which at 31 deg is **13 315 s =
66 500 iterations** at `dt` = 0.2 s. This ran 100, i.e. **0.15 % of one inertial time unit**, and
showing a circulation response would take ~155 hours of wall clock at the current step.
**So the term is verified CONNECTED, CORRECTLY SIZED and STABLE, and it is NOT yet verified to do
anything to the climate.** The obvious lever is `dt`, exactly as this file concluded for
`ATM_SFC_FLUX` — one hour of physical time is 17 973 iterations. **THAT LEVER HAS NOW BEEN TRIED
AND IT DOES NOT REACH: the next subsection measures the `dt` ceiling at ~4x on the model's own
precipitation, where `1/f` needs 16x and 16 600 iterations.**

**DEFAULT STAYS 0.0** for that reason and no other: nothing about it has been shown on a spun-up
field, and this tree flips defaults on measurements, not on arguments. **The attempt to reach the
spun-up measurement by raising `dt` failed, and failed in an instructive way — see the next
subsection before repeating it.**
### `dt` CANNOT BE RAISED FAR ENOUGH TO REACH THE ADJUSTMENT, AND THE LADDER THAT SAID IT COULD WAS JUDGED ON THE WRONG CRITERION

**THE ATTEMPT.** Geostrophic adjustment takes `1/f` = 13 315 s at 31 deg = **66 500 iterations** at
`dt_visc` = 1e-4, so `ATM_HYDRO_PGF`'s effect on the CIRCULATION is unreachable at the shipped
step. `dt_visc` was swept 2x, 4x, 8x, 16x, 60 iterations each from the iteration-600 checkpoint,
and **all four passed exit 0, zero NaN, and `max|u|` within 5 % of the control.** On that basis a
paired 1000-iteration run (`ATM_HYDRO_PGF` on / off, `dt_visc` = 16e-4 = 3.20 s/iter) was launched
and reached iteration 1600 clean.

**THE LADDER WAS WRONG, AND THE MODEL'S OWN PRIMARY DIAGNOSTIC SAYS SO.** Precipitation at
iteration 660, which nobody looked at until the long run had finished:

| `dt_visc` | s/iter | Precip mm/a | vs 1e-4 | pattern r | `max|w|` |
|---|---|---|---|---|---|
| **1e-4 (shipped)** | 0.20 | **1014.9** | — | +0.458 | 17.97 |
| 2e-4 | 0.40 | 1019.4 | +0.4 % | +0.454 | 18.27 |
| 4e-4 | 0.80 | 1048.4 | +3.3 % | +0.456 | 18.57 |
| 8e-4 | 1.60 | 1127.3 | **+11 %** | +0.454 | 19.17 |
| **16e-4** | 3.20 | **1598.1** | **+57 %** | +0.422 | 20.35 |

And at iteration 1600 the 16e-4 pair reports **8733 mm/a — 793 % of NASA** — with `r` degraded
0.458 -> 0.263 and `max|w|` **90.2 m/s** in the hydrostatic arm against a +-100 clamp. **The runs
never NaN'd and were never valid.**

**SO THE USABLE CEILING IS ~4x, NOT 16x**, and that does not help: `1/f` is 16 600 iterations at
`dt` = 4e-4, ~46 hours of wall clock, and at 4x the per-iteration filters are already 4x weaker per
unit physical time. **"Raise `dt` and run to adjustment" is not reachable in this model.**

*This is the tree's own recorded trap, and it was walked into with the warning in the file.*
CLAUDE.md already says **"Check `Precip mean` against the NASA line before calling any moisture
change an improvement"**, recorded as the sixth occurrence of "the answer was in output nobody
opened". A `dt` change is not a moisture change, which is exactly why the rule was not applied —
and `dt` multiplies every rate in the microphysics. **Exit 0 and no NaN is not a stability
criterion in this tree.** Seventh occurrence.

**TWO CLAIMS MADE FROM THAT RUN ARE RETRACTED.** An interim reading at iteration 660 reported
`Psi(ground)` **-14.9 %** with its growth arrested. Neither survives: by iteration 1200 the sign
reverses and at 1600 the hydrostatic arm is **+22.6 %** ABOVE the control (4.734e+11 against
3.863e+11), and the run it came from is invalid anyway. **`ATM_HYDRO_PGF`'s effect on
`Psi(ground)` is UNMEASURED.**

**AND THE ONE THING THE INVALID RUN SUGGESTED IS TRUE FOR A REASON THAT NEEDS NO RUN AT ALL.**
In the hydrostatic arm `vbar` at the jet core did not settle toward balance — it crossed zero near
iteration 1100 and ran away the other way, -0.21 -> -0.20 -> **+3.40**, while the control ran away
negatively, -0.21 -> **-3.46**. That is what the structure requires: **`p_stat` is DIAGNOSED from
the temperature** — `densities()` builds it as an analytic barometric function of layer height and
`t.x[0]` — **and it does not respond to the flow's mass convergence at all.** So a prescribed
hydrostatic pressure gradient is a FORCING, not an adjustment mechanism: there is no restoring
feedback, and `v` integrates whatever residual imbalance is left rather than settling into
geostrophy. Real adjustment needs the pressure to answer the flow, i.e. a prognostic surface
pressure or a mass-continuity equation, and this model's only flow-responsive pressure is `p_dyn`
— which the section above measures at four to five orders too weak.

**WHAT THAT MEANS FOR THE KNOB.** `ATM_HYDRO_PGF` can hold the model NEAR thermal-wind balance
where the temperature field is right; it cannot drive the model TOWARD it. The measurements that
stand are the ones taken at the shipped `dt` in the subsection above — connected, correctly sized,
stable, byte-identical off, and the first geostrophically balanced cells this model has had. The
circulation response remains unmeasured, and the obstacle is no longer "run it longer": it is that
this model has no pressure that answers to the flow.
### The balanced initial state DOES give this model a mid-latitude jet — and the precipitation cannot feel it, for a reason that retires the whole class of experiment

**`ATM_TW_BALANCE=1.0` AGAINST A CLEAN CONTROL, 600 ITERATIONS FROM SCRATCH EACH, `config_accept.xml`
DEFAULTS, MOIST PHYSICS FROM ITERATION 0, 24 THREADS, ONE BINARY, ~61 MIN PER ARM** (2026-09-03,
`output_twctl` / `output_tw10`; both exit 0 with **zero NaN**; the two `[RUN CONFIG]` banners differ
in `TW_BALANCE` and nothing else — `HYDRO_PGF=0` in BOTH, so this is the initial state alone).

From scratch is mandatory rather than preferred: `load_state()` overwrites everything the init path
produces, so a restart arm cannot exercise this knob at all — and the byte-identical check that
cleared it was a restart, i.e. it could only ever have verified the OFF branch.

**THE CONTROL IS ALSO THE FROM-SCRATCH OFF-BRANCH CHECK, AND IT REPRODUCES `output_rsfull100` TO
EVERY RECORDED DIGIT** — zonal-mean jet above 3 km **26.4877 / 22.5733 / 20.4221 / 18.3027** at
iterations 20 / 160 / 320 / 600 against the recorded 26.488 / 22.573 / 20.422 / 18.303, core at
5513 m. That binary has since gained `ATM_TW_BALANCE`, `ATM_EVAP_SPREAD`, the pressure-clamp
instrument and the re-sized clamps; none of them touches the default branch.

**THE KNOB DOES EXACTLY WHAT IT WAS WRITTEN TO DO.** At initialisation: 64 619 columns, **largest
shear added 242.45 m/s**, capped at the 80 m/s ceiling in 6798 cells. Zonal-mean jet above 3 km:

| iteration | 20 | 200 | 400 | **600** | core at 600 | latitude |
|---|---|---|---|---|---|---|
| **control** | 26.488 | 21.919 | 19.651 | **18.303** | 5513 m | **31S** |
| **`TW_BALANCE=1.0`** | 45.261 | 37.130 | 33.428 | **31.232** | 5513 m | **65S** |

**A 71 % STRONGER JET, AND IT IS IN THE MID-LATITUDES WHERE THE CONTROL HAS NOTHING.** The control's
zonal-mean maximum sits at 30-31S for the whole run; this one sits at 64-65S from iteration 20
onward. That is a storm-track jet, in the band the model is 6x too dry in.

**AND THE DECAY IS THE SAME FRACTION, WHICH IS THE THIRD CONFIRMATION THAT THE FILTER IS LINEAR.**
-30.9 % for the control against **-31.0 %** for the balanced arm over the same 600 iterations, and
`dw_radial` in the 20-70 deg band above 3 km is 2.09e-03 against **4.38e-03**, 2.1x — tracking the
amplitude, as a linear damping must. The balanced state is not held by anything: with `HYDRO_PGF=0`
there is no pressure gradient to balance it, so it is a larger initial condition being eroded at the
same rate, exactly as the strength sweep predicted.

**THE EDDIES DO RESPOND. THE PRECIPITATION DOES NOT.** Band-median eddy momentum flux convergence
`wbud_advh` goes 6.5e-08 -> **2.15e-07, 3.3x** — still 20x below Coriolis, but plainly alive. Then,
from the two runs' own written fields at full precision (not the 1-decimal print):

| cos-lat mean, mm/a | control | `TW_BALANCE=1.0` | change | NASA |
|---|---|---|---|---|
| 0-15 | 3436.381 | 3440.865 | +0.131 % | 1487.0 |
| 15-35 | 273.961 | 273.142 | -0.299 % | 761.4 |
| **35-65** | **156.628** | **156.649** | **+0.013 %** | **981.1** |
| 65-90 | 7.481 | 7.474 | -0.088 % | 364.2 |
| precipitable water | — | — | **+0.000 % in every band** | |
| surface temperature | — | — | **within 0.015 %** | |

**156.6 FOR THE FOURTH TIME, AND IT IS NOT A PRINT-PRECISION ARTEFACT** — 0.013 % is measured off the
fields. The Shapiro sweep got 156.6 at two filter strengths; this gets 156.6 with the jet 71 %
stronger AND relocated 34 degrees poleward into the band itself.

**THE REASON IS ARITHMETIC, IT WAS AVAILABLE WITHOUT RUNNING ANYTHING, AND IT APPLIES TO EVERY
DYNAMICAL ARM IN THIS FILE.** 600 iterations is **120.2 s** of physical time. The 35-65 deg band's
median \|`vbar`\| above 500 m is **0.33 m/s**, so a parcel moves **36 metres — 0.03 % of a 1-degree
cell.** One cell meridionally is **1.85 MILLION iterations**; even the 31 m/s jet needs **18 500**
to move air one cell zonally. **Precipitation here is a local column response to a thermodynamic
state that transport cannot have touched** — which the table above measures directly: precipitable
water identical to six figures in every band, temperature to 0.015 %, and the precipitation moving
by about as much as the local moisture did and no more.

**SO THE FOUR "THE JET BUYS NOTHING" ARMS ARE ONE MEASUREMENT, AND IT IS OF THE ADVECTIVE
DISPLACEMENT, NOT OF THE EDDIES.** *"The eddies are STARVED, not damped"* stands where it was
measured — on `wbud_advh` against `wbud_cor` in the momentum budget — and **must not be read off the
precipitation bands, which cannot respond to either.** Likewise *"do not expect a filter setting to
fix the 35-65 degree band"* generalises: **do not expect ANY dynamical change to move ANY
precipitation band in a run this tree can afford.** A circulation experiment scored on precipitation
is scoring nothing, whatever the size of the circulation change.

**WHAT IS AND IS NOT SETTLED.** Settled: the balanced initial state is written, stable at full
strength, byte-identical off, from-scratch verified, and it is the only thing in this tree that has
ever put a jet in the mid-latitudes. Not settled: whether that jet does anything to the climate,
because no run of this length can show it. The lever is the same one the `dt` subsection above
failed to reach — physical integration length — and the honest next instrument is a
circulation-side score (eddy heat/momentum flux, `Psi`, the jet itself) rather than another
precipitation table. `Psi(ground)` is -0.4 % between the arms, as it must be: `ATM_TW_BALANCE_V` is
off, so only the ZONAL wind was touched and `Psi` is meridional.

### What this does and does not say

It IS: the meridional pressure gradient measured against Coriolis on a spun-up default-
configuration field; the mechanism read off the source; the leading candidate repair swept and
retired; the pressure's normalisation corrected; and the first of the three repairs written and
measured (the subsection above). **It is NOT yet a repair of the CLIMATE**: `ATM_HYDRO_PGF` is
verified connected, correctly sized, stable at full strength and byte-identical off, and it moves
no velocity on any run this tree can afford. **ALL THREE ROUTES HAVE NOW BEEN WRITTEN AND
MEASURED, AND NOT ONE MOVES A PRECIPITATION BAND.** `ATM_HYDRO_PGF` (force side, jet-core
ageostrophic residual 0.9999 -> 0.2498); `ATM_BUOY_CONSISTENT` (the buoyancy's missing factor of
5.0e5, reaching momentum through the model's own elliptic pressure, band p05 residual 0.996 ->
0.437 at 100 iterations and 0.336 at 600 — and it did NOT turn out to be unstable, as this file
expected); and `ATM_TW_BALANCE`, the balanced initial state in the manner of ATHAD's
`initBalancedState`, which is the only one that moves a VELOCITY — a 71 % stronger jet, relocated
into the storm-track band — and moves the 35-65 deg precipitation by **0.013 %**. The subsection
above says why, and the reason is not a property of any of the three.
And one more caveat with teeth: **one iteration is 0.2 s, so 1/f at 31 deg is 13 315 s = 66 500
iterations.** The longest run in this tree is 1000. Nothing here has been integrated for 2 % of an
inertial time unit — which does NOT excuse the missing balance (a projection pressure is
diagnostic and re-solved every step, so it does not need an inertial period to appear), but does
mean no arm in this file can be read as an adjusted state.

## The ocean has the atmosphere's defect one model over, and `HYD_PHYDRO_SALT` repairs the field the missing force would need

**`RHS_Hyd_Turb.cpp:722` HAS BEEN SAYING THIS IN A COMMENT ALL ALONG**: *"a real baroclinic ocean
needs the horizontal grad(p_hydro) in rhs_v/rhs_w, not just this term."* The ocean's only route
from temperature or salinity into momentum is `buoy_nd` in `rhs_u` — RADIAL only, and the same
comment records it as carrying "the spurious dt of the pre-`5f71571` convention (~1e-7 too small),
so buoyancy is effectively OFF and the ocean is dynamically barotropic (wind + Coriolis only)".
**That is `ATM_HYDRO_PGF`'s finding in the other model**: no horizontal pressure gradient in the
horizontal momentum equations, and the one vertical surrogate inert.

**AND THE FIELD THAT FORCE WOULD BE BUILT FROM DOES NOT KNOW ABOUT SALT.** `p_hydro` integrates
`r_water` — thermal expansion and compressibility — while `r_salt_water`, computed three lines
below in the same loop, is the full Gill approximation in temperature, SALINITY and depth.
`HYD_PHYDRO_SALT=<0|1>` (**default 0 = shipped**, first environment knob in the hydrosphere)
integrates the latter instead, with a plausibility floor and an unconditional one-line print.

**FOUR ARMS. THE TWO THAT MATTER ARE AT ONE THREAD, BECAUSE THIS MODEL IS NOT BIT-REPRODUCIBLE
UNDER OpenMP** (`project_hydro_nondeterminism`: a race survives in the momentum path). A
multi-threaded pair cannot tell "the knob changed nothing" from "the threads did".

| arm | binary | threads | nm |
|---|---|---|---|
| `qb0` | `hyd_base`, HEAD, knob absent | 1 | 20 |
| `qn0` | new binary, knob unset | 1 | 20 |
| `qn1` | new binary, `HYD_PHYDRO_SALT=1` | 1 | 20 |
| `ps0` / `ps1` | new binary, off / on | 24 | 300 |

All four exit 0 with **zero NaN**, on the current default atmosphere's own surface forcing
(`output_twctl/0Ma_smooth_Transfer_Atm_600.vwtp`, yesterday's from-scratch control).

**OFF-BRANCH: BYTE-IDENTICAL OVER ALL 19 WRITTEN FILES**, `qb0` against `qn0` — every VTK slice,
the panorama `.vts`, both SST transfer files, `PlotData_Hyd.xyz`, all four budget CSVs, and the
**364 MB restart binary** (`md5 f8de6005...` both). Not one byte.

**AND THE KNOB IS DYNAMICALLY INERT, MEASURED RATHER THAN READ OFF THE SOURCE.** `qn1` against
`qn0` at one thread: **16 of the 17 restart arrays are BIT-IDENTICAL — `t`, `u`, `v`, `w`, `c`,
their `n` copies, `p_dyn`, `tke`, `dis`, `nue` all differ by exactly 0.0** — and only `p_hydro`
moves, by 2.30 % rms. That is the dependency graph confirmed: `p_hydro` and `r_water` reach only
`PresGradForce` and `CentrifugalForce`, which are diagnostics **nothing in `rhs_u`/`rhs_v`/`rhs_w`
reads**, plus the restart and the VTK writers. `r_salt_water` is not touched at all — its depth
term is `p_km = d_i*step_km`, GEOMETRIC depth, not `p_hydro`.

**So this knob cannot change the ocean today, and that is the point of it rather than a
disappointment**: it repairs the field, and the force that would consume it does not exist yet.

**THE MAGNITUDE, ON THE 300-ITERATION 24-THREAD PAIR** (grid 41x181x361, and note the ocean is in
SHALLOW mode: `p_hydro` maxes at 20.6 bar, a ~200 m column). Cos-lat mean `p_hydro`
**10.9125 -> 11.1651 bar**; over fluid cells the rms difference is **2.29 %**, max 0.60 bar. Every
other field in the restart moves 1e-4 to 1e-9 relative — the documented non-determinism, worst in
the velocity path exactly where that memory says the surviving race is.

**BUT THE MEAN IS NOT THE QUANTITY. A CONSTANT OFFSET IS INVISIBLE TO A GRADIENT, AND THE GRADIENT
IS THE FORCE.** Median \|(1/rho) dp/dy\| over the extratropical fluid interior:

| | median \|(1/rho)dp/dy\| | \|pgf\|/\|f w\| median | rms change |
|---|---|---|---|
| **all fluid cells** shipped | 5.33e-08 | 0.126 | — |
| all fluid cells, salt-aware | **8.42e-07** | **2.105** | 2.0 % of the shipped rms |
| **salty columns only** shipped | 5.71e-08 | 0.157 | — |
| salty columns, salt-aware | **6.27e-07** | **1.619** | 3.3 % of the shipped rms |

**AN ORDER OF MAGNITUDE, NOT A CORRECTION — AND IT IS NOT THE SALINITY BUG.** Restricting to
columns carrying no fresh cell at any level above (63.8 % of the scored ocean) it is still
**11.0x**. The rms change is only 2-3 % because the rms is set by a few coastal cells where both
branches are large; the TYPICAL cell is where the shipped field has essentially no horizontal
structure at all.

**AND THE CHANGE'S OWN RATIONALE IS HALF WRONG: IT IS NOT MAINLY THE SALT.** The knob's comment
argues from `beta_p` ~ 0.8 kg/m3 per psu against `alfa_t_p` ~ 0.2 per degC. Rebuilding the column
integral from the Gill equation with one variable frozen at its global fluid mean at a time
(post-processing of `ps0`'s own restart, no re-run), median \|dp/dy\| in salty columns:

| | thermal only (S = S_bar) | haline only (T = T_bar) | both | neither |
|---|---|---|---|---|
| median \|dp/dy\| | **6.83e-07** | 3.71e-07 | 6.20e-07 | **0.0** |
| lateral sigma(rho) | 2.22 | 9.92 | 8.68 kg/m3 | 0.00 |

**The THERMAL half is the larger one**, the two partly oppose (both < thermal alone), and the
`neither` row returning exactly 0.0 is the pipeline's own control. **The reason is that shipped
`r_water` has no absolute temperature in it either**: it is anchored at
`r_water[im-1] = r_0_water`, a GLOBAL CONSTANT, and propagated downward by *increments* —
`beta_water*(T_i - T_{i+1})` between adjacent LEVELS, plus compressibility. A density built from
vertical differences of a constant carries no lateral structure whatever the temperature does. So
the shipped `p_hydro` is **barotropic by construction**, and what the knob restores is the whole
equation of state, of which salinity is the smaller half.

**AND THE PLAUSIBILITY FLOOR IS CATCHING 12.35 % OF THE OCEAN, NOT THE ~3 % ITS COMMENT CLAIMS.**
Measured on `ps0`'s own field: **11.62 % of fluid cells are below 5 psu and 28.75 % below 30 psu**,
the floor rejects **12.35 %** (the model's own print says 214 661 per call, against 222 850
reconstructed — 3.7 % apart, which is the reconstruction's flat-depth approximation), and the
rejected cells have median salinity **0.028 psu** and span **every level, i = 0 to 40** — not the
coast-surface artefact the recorded bug describes. `project_hydro_salinity_ic`'s "coast-adjacent
water cells still 0" is a much larger defect than that phrase suggests, and it is now the largest
single obstacle to using this field.

**TWO LIMITS OF THE FLOOR, BOTH MEASURED RATHER THAN ARGUED.** It is a fixed 1005 kg/m3 threshold,
and `C_p` alone reaches 1019 at 4 km, so on a deep grid it would stop catching fresh cells
entirely; it works here only because the column is 200 m. And falling back to `r_water` at a
rejected cell puts a density DISCONTINUITY next to its salty neighbour, which manufactures exactly
the horizontal gradient this field is being repaired to provide — which is why the all-cells 15.8x
must be read as the salty-column **11.0x** plus contamination, and why the salinity IC is a
prerequisite rather than a follow-up.

**WHAT THIS IS AND IS NOT.** It IS: the field a horizontal baroclinic PGF would be built from,
repaired, verified inert, and measured at an order of magnitude on the quantity that force would
use. It is **NOT a change to the ocean's circulation and cannot be one** — `p_hydro` and `r_water`
reach only `PresGradForce` and `CentrifugalForce`, which are diagnostics nothing reads, plus the
restart and the VTK. **The force itself is still unwritten.** Default stays 0 for the same reason
`ATM_HYDRO_PGF`'s does: nothing about it has been shown on a circulation, and this tree flips
defaults on measurements.

## Open risks

- **`ATM_CLOUD_FRAC`: the sub-grid cloud scheme is WRITTEN AND STRUCTURALLY RIGHT, AND IT IS NOT
  CALIBRATED. DO NOT FLIP IT.** (2026-08-29, default 0, off-branch unchanged.)

  Uniform-PDF (Smith / Sundqvist) closure, one parameter, the `H_crit` already present. Total
  water uniform over `[qbar - D, qbar + D]` with `D = (1 - H_crit)*q_sat`, cloud where
  `q_t > q_sat`:

      s   = qbar + D - q_sat = q_sat*(RH - H_crit)
      f   = clamp(s / 2D, 0, 1)
      q_c = f^2 * D        (f < 1),      qbar - q_sat   (f = 1)

  `f` = 0 exactly at `RH = H_crit`, `q_c` is the GRID-MEAN condensate so the in-cloud value
  `q_c/f` stays finite as `f -> 0`. Computed on the fly, no new field, so
  `sizeof(cAtmosphereModel)` is untouched and the stack-canary hazard does not apply.

  **WHAT WORKS.** On the shipped humidity it does exactly what it should to the STRUCTURE:
  cloudy cover **98.94 % -> 22.16 %** and the column path 1584 -> 4.59 g/m2. The pathology of
  "cloud in 99 % of columns over 38 of 41 levels" is gone.

  **AND `ATM_RH_PROFILE` GIVES A REALISTIC HUMIDITY**, which is the other half of the repair:

  | z | 0 m | 441 | 1368 | 3316 | 4988 |
  |---|---|---|---|---|---|
  | mean RH, shipped | 0.933 | 0.813 | 0.896 | 0.929 | 0.937 |
  | mean RH, `ATM_RH_PROFILE=1` | **0.719** | **0.610** | **0.601** | **0.485** | **0.392** |

  0.72 at the surface falling to 0.39 at 5 km, against Earth's ~0.8 falling to ~0.4-0.5. **That
  is the first physically-shaped humidity profile this model has had.**

  **WHAT DOES NOT WORK, AND IT IS UNRESOLVED.** Put the two together and the condensate
  vanishes rather than landing at the observed 50-100 g/m2:

  | `ATM_RH_CRIT` (with `RH_PROFILE` + `CLOUD_FRAC`) | 0.55 | 0.40 | 0.30 | 0.20 |
  |---|---|---|---|---|
  | column path g/m2 | 0.0006 | 0.0043 | 0.151 | **0.627** |

  **THE ~100x DISCREPANCY IS FOUND, AND IT IS A GRID-MEAN SATURATION ADJUSTMENT DESTROYING A
  SUB-GRID CLOUD.** The arithmetic was right and the sweep was measuring an already-destroyed
  field. Instrumented at both ends (`ATM_CLOUD_INIT_DIAG` now prints the same column path
  `ATM_CWP_CENSUS` reports, at init and again after the init adjustment):

  | | column path at init | after the init `SaturationAdjustment` | removed |
  |---|---|---|---|
  | shipped (near-saturated column) | 2547.42 | 1588.65 | 38 % |
  | `RH_PROFILE` + `CLOUD_FRAC` | **270.70** | **0.633** | **99.77 %** |

  **`initCloudIce` with the fractional closure produces 270.70 g/m2 — a physically sensible
  number, ~3x the observed 50-100 rather than 30x.** The scheme works. What follows it does not:
  `cAtmosphereModel.cpp:461` calls `SaturationAdjustment(*this).run()` **UNCONDITIONALLY at
  initialisation**, outside the `moist_phys_start_iter` gate, and it drives `q_v` toward
  `q_sat` on the GRID MEAN. A fractional scheme puts condensate exactly where the grid mean is
  SUBSATURATED, so a grid-mean adjustment evaporates precisely the cloud the closure created —
  a factor of **428**.

  **AN EARLIER TEST OF THIS WAS INVALID AND SAID SO WRONGLY.** Setting the moist gate to 300
  changed nothing (0.6266 vs 0.6288) and that was read as exonerating `SaturationAdjustment`.
  The gate governs only the IN-LOOP call; the init call at line 461 is ungated and is the one
  that does the damage.

  **SO THERE ARE THREE MUTUALLY-CONSISTENT WRONG CHOICES, NOT TWO**: a constant near-saturated
  RH, an `H_crit` of 0.80-0.98 that only such a column can exceed, and a grid-mean saturation
  adjustment that only such a column survives. The shipped model is self-consistent and every
  one of the three is wrong. Fix any one alone and the condensate collapses; the 38 % vs 99.77 %
  row above is that statement measured.

  **WRITTEN, AND THE CHAIN NOW CLOSES.** `SaturationAdjustment` targets the FRACTIONAL
  equilibrium under `ATM_CLOUD_FRAC`: with total water `q_t = q_v + q_c + q_i` — conserved by
  the loop, since `d_cnd + d_dep = d_q_v` — the target becomes `q_t - q_c_eq` with `q_c_eq` from
  the same uniform-PDF closure. It reduces EXACTLY to the shipped `q_s` when `H_crit` = 1
  (`D` = 0, `f` = 1), so the off-branch is unchanged BY CONSTRUCTION rather than by a guard, and
  measures so (1583.9636 against a 1583.9636 baseline). The entry test is widened too: a cell
  can now be cloudy while the grid mean is subsaturated, so admission cannot be conditioned on
  `q_v > q_sat` alone.

  **THE CONDENSATE NOW SURVIVES.** `nm = 20`, `ATM_RH_PROFILE` + `ATM_CLOUD_FRAC`:

  | `ATM_RH_CRIT` | at init | after init SatAdj | at the radiation | columns w/ cloud |
  |---|---|---|---|---|
  | 0.35 | 26.81 | 25.09 | **24.99** | 73 % |
  | 0.32 | 55.34 | 53.61 | **53.39** | 82 % |
  | 0.30 | 80.26 | 79.34 | **79.02** | 86 % |
  | 0.28 | 109.77 | 110.59 | **110.14** | 89 % |

  **Loss across the adjustment is 1.5 %, against 99.77 % before** — that single row is the
  repair. And the column path is now a CALIBRATABLE quantity that lands in the observed
  50-100 g/m2 at `ATM_RH_CRIT` ~ 0.30-0.32, instead of being forced there by `cwp_cap_col`
  dividing by 79 three modules downstream.

  **TWO CAVEATS ON READING THAT TABLE.** The last column is the fraction of COLUMNS carrying
  more than 5 g/m2, NOT satellite cloud cover — a column with `f` = 0.1 at a few levels counts
  as cloudy — so 86 % is not a claim of 86 % cloud cover and must not be compared with Earth's
  ~67 %. And `nm` = 20 is the INITIAL field: this is the condensate the model starts from, not
  one a spun-up circulation maintained.

  **ACCEPTED AND FLIPPED ON BY DEFAULT 2026-08-31.** Eleven defaults moved together, because
  every one of them is wrong alone: `ATM_RH_PROFILE` 0 -> 1, `ATM_RH_CRIT` 0.8 -> 0.30,
  `ATM_CLOUD_FRAC` 0 -> 1, `ATM_CWP_CAP` 20 -> DISABLED, `ATM_CLOUD_RAD_FRAC` 0 -> 1,
  `ATM_QC_CRIT` 0.5 -> 0.05 g/kg, `ATM_ICE_COLD` 0 -> 1, `ATM_T_FLOOR` 236.15 -> 216.65 K,
  `ATM_RH_MIN` 0 -> 0.65, `ATM_RH_MIN_LAT` 0 -> 1, `ATM_RH_MIN_PTOP` 0 -> 475 hPa.
  **Setting each variable back to the value above restores the old branch**, and both directions
  are verified at 24 threads: a clean environment reproduces the `c475` arm (OLR clear 263.404
  vs 263.396, cloudy 233.228 vs 233.227, temperature extremes bit-identical) and the full revert
  reproduces the shipped branch (precip 1048.56 vs 1048.54, OLR 273.08818164 vs 273.08818780,
  precipitable water 49.64 both, min T -37.151524 vs -37.151522) -- the residuals are the
  documented fixed-thread-count non-determinism, not a configuration difference.
  **TWO THINGS THE ACCEPTANCE DOES NOT CLAIM.** It is measured at `nm` = 100, which is 20
  SECONDS of physical time, so it is the initial field plus a short transient and not a spun-up
  climate. And the tropics are 100 % cirrus-covered against Earth's ~40 %, which comes WITH the
  configuration rather than being something left to tune: a zonally uniform floor cannot make
  longitudinally patchy cirrus.

  **AND THE FIRST LONG RUN SAYS THE PRECIPITATION IS NOT CONVERGED** (2026-08-31, restart from
  the accepted configuration's own iteration-100 state, run to 1000 at 24 threads, 41 + 51 min;
  no NaN anywhere, `max|w|` 18.5 -> 16.8 m/s against the +-100 clamp, so the old MoistConvection
  runaway does NOT return). Precipitation, high parity of the 2dt sawtooth:

  | iteration | 100 | 200 | 400 | 600 | 700 | 800 | 900 | 1000 |
  |---|---|---|---|---|---|---|---|---|
  | Precip mm/a | 992 | 988 | 979 | 1003 | 1021 | 1087 | 1141 | **1215** |
  | low parity | 335 | 331 | 337 | 331 | 374 | 429 | 546 | **713** |

  **FLAT TO ITERATION ~600, THEN RISING: +22 % on the high parity and +113 % on the low one over
  the last 400 iterations, and still climbing at the end.** The sawtooth is collapsing as the low
  parity catches up, which is what a growing signal does to a 2dt alternation.

  **THE CLOUD AND THE RADIATION ARE STABLE WHILE THE RAIN IS NOT**, and that is the shape of the
  problem: LWP 103.0 -> 107.5 -> 116.3 (+13 %), IWP 20.09 -> 20.30 -> 20.65 (+3 %), RH at 10.9 km
  44.4 -> 44.4 -> 44.5 %, precipitable water 30.27 mm unchanged, clear OLR 263.40 -> 263.33,
  all-sky 233.23 -> 233.07 at iterations 100 / 600 / 1000. **The precipitation is growing four
  times faster than the condensate that feeds it**, so this is a conversion-efficiency drift, not
  a moistening. `P_conv` also collapses over the same stretch, 0.32 -> **0.0196 mm/a**, while
  `P_rain` goes 970 -> 1197: the stratiform path is taking over from the convective one.

  **SO THE 992 mm/a THAT THE CONFIGURATION WAS ACCEPTED ON IS AN ITERATION-100 NUMBER, AND
  ITERATION 100 IS 20 SECONDS.** At 1000 iterations -- 200 seconds -- it is 1215, 24 % above
  NASA. Nothing here is spun up; the caveat recorded at acceptance is now measured rather than
  merely stated. **Do not quote the accepted configuration's precipitation without its iteration
  count.**

  **The paragraph below is the pre-acceptance record, kept because it is what the flip was
  judged against.** STILL DEFAULT OFF, ALL FOUR (`ATM_RH_PROFILE`, `ATM_RH_CRIT`, `ATM_CLOUD_FRAC`, and the
  fraction-aware adjustment they share). What has NOT been done: `cwp_cap_col` is still 20 and
  still divides the now-physical path by ~4, so it must be disabled in the same flip; the
  radiation still treats every column as overcast rather than weighting by `f`; and none of this
  has been run to a spun-up state or checked against the OLR.

  **So the state is: the humidity fix is measured and good, the closure is measured and
  structurally good, and the two together are NOT yet a working pair.** All three knobs
  (`ATM_RH_PROFILE`, `ATM_RH_CRIT`, `ATM_CLOUD_FRAC`) stay default OFF. The next step is the
  ~100x discrepancy above, not another sweep and not a flip.

  **THE FULL FLIP HAS NOW BEEN RUN TO `nm` = 400 WITH MOIST PHYSICS FROM ITERATION 0, AND IT
  DESTROYS THE PRECIPITATION. DO NOT FLIP IT** (2026-08-30, 24 threads,
  `ATM_RH_PROFILE=1 ATM_RH_CRIT=0.30 ATM_CLOUD_FRAC=1 ATM_CWP_CAP=off ATM_CLOUD_RAD_FRAC=1`,
  against the identical config shipped).

  | | shipped | full flip | Earth / NASA |
  |---|---|---|---|
  | column condensate path | 1583.2 g/m2 | **78.8** | 50-100 |
  | p05 / p50 / max | 213 / 1959 / 2298 | **1.0 / 86 / 200** | — |
  | levels carrying cloud | 38 of 41 | 18 of 41 | — |
  | precipitable water | 49.6 mm | **29.6** | ~25 |
  | **Precip mean** | **1046 mm/a** | **1.21 mm/a** | **978 (NASA)** |
  | P_rain mean | 805 | **0.0092** | — |
  | P_snow mean | 228.6 | 1.22 | — |
  | P_conv mean | 9.65 | **0.000** | — |
  | cloud LW forcing | 25.35 W/m2 | **9.05** | ~25 |
  | all-sky OLR | 247.7 | 266.5 | ~240 |

  **THE CONDENSATE AND THE HUMIDITY BOTH BECOME RIGHT AND THE MODEL STOPS RAINING.** 78.8 g/m2
  is in the observed band and 29.6 mm of precipitable water is better than the shipped 49.6, on a
  SPUN-UP field at iteration 400 -- the "`nm` = 20 is the INITIAL field" caveat above is
  discharged. And precipitation falls by **865x**, rain by **87 500x**, convective precipitation
  to exactly zero.

  **THE CAUSE IS A HARD THRESHOLD FITTED TO THE BROKEN CONDENSATE**, `TwoCatIceScheme.h:294`:

      constexpr double q_c_crit = 5.0e-4;            // 0.5 g/kg Kessler threshold
      if(m.cloud.x[i][j][k] > q_c_crit)
          S_c_au = c_c_au * (m.cloud.x[i][j][k] - q_c_crit);
      else S_c_au = 0.0;

  Under the flip the grid-mean cloud water at the profile peak is ~0.04 g/kg -- **twelve times
  below the threshold** -- so `S_c_au` is identically zero everywhere and no rain forms at all.
  Snow survives only because ice autoconversion has no threshold (`if(ice > 0.0)`). The
  threshold's own comment records that it was installed because without it "cloud~0 yet huge
  P_rain -> gross precip ~30x NASA": **it was fitted against the 20x-too-large condensate.**

  **SO THE SHIPPED 1046 mm/a IS A FIFTH CANCELLING PAIR** -- an over-large condensate feeding a
  threshold raised to compensate -- and it is the one that matters most, because it is the
  model's headline output and it currently matches NASA to 7 %.

  **TWO THINGS ARE WRONG WITH `q_c_crit` UNDER A FRACTIONAL SCHEME, AND ONLY THE SECOND IS A
  TUNING QUESTION.** (1) It is applied to the GRID MEAN. Autoconversion is an in-cloud
  microphysical process and must see `q_c/f`. **This is structurally the same defect as the
  grid-mean `SaturationAdjustment` fixed in `6d163a0`, one module downstream** -- the fractional
  correction was made in the adjustment and not in the microphysics that follows it. (2) 0.5 g/kg
  is defensible in-cloud and never occurs as a grid mean over a 1x1 degree box.

  **SO THE FLIP LIST IS INCOMPLETE, NOT MERELY UNCALIBRATED.** A fraction-aware autoconversion in
  the ice scheme is a PREREQUISITE, not a follow-up. Until it exists the flip cannot be judged:
  its condensate, humidity and cloud-cover fields are all better than the shipped ones and its
  precipitation is three orders of magnitude wrong.

  **AND THE CLOUD LW FORCING FALLS 25.35 -> 9.05 W/m2 FOR A SECOND, INDEPENDENT REASON: THE FLIP
  REMOVES THE HIGH CLOUD.** Per-level condensate, cos-lat mean g/m2:

  | z | 2163 m | 3316 m | 4988 m | 7412 m | 10927 m |
  |---|---|---|---|---|---|
  | shipped | 99.8 | 127.2 | 92.8 | 30.8 | **20.5** |
  | flip | 2.4 | 10.6 | 6.7 | 0.03 | **0.000** |

  Earth's ~25 W/m2 of LW cloud forcing is mostly cold high cloud; this leaves only low cloud,
  which radiates near surface temperature and forces almost nothing. `ATM_RH_PROFILE`'s
  Manabe-Wetherald profile falls too fast aloft for cirrus to form. That is a THIRD prerequisite,
  independent of the autoconversion one.

  **THE AUTOCONVERSION IS NOW FRACTION-AWARE (`3ea78e3`), IT BOUGHT x1850 ON RAIN, AND THE
  THRESHOLD IS NOW EXHAUSTED AS A LEVER** (2026-08-30, `nm` = 100, gate 0, 24 threads).
  Autoconversion is a LOCAL process, so the grid-mean tendency is `f*R(q_c/f)`. Only the
  NONLINEAR terms move under that transform -- `f*c*(q_c/f) = c*q_c` identically -- and every one
  that moved is nonlinear because it carries a threshold: TwoCat `S_c_au`/`S_ac`, ZeroCat `S_au`,
  OneCat `S_au`/`S_rim`/`S_shed`. `S_rim`, `S_shed`, `S_c_frz`, `S_i_au`, `S_nuc` are linear and
  correct as they stand.

  The structural repair alone takes rain **0.0092 -> 17.05 mm/a**. Then `ATM_QC_CRIT` (new,
  g/kg, default 0.5 = shipped):

  | `q_c_crit` g/kg | 0.20 | 0.15 | 0.10 | 0.05 | 0.02 | 0.01 | NASA |
  |---|---|---|---|---|---|---|---|
  | Precip mm/a | 140/370 | 159/441 | 181/510 | 207/563 | 214/618 | **222/638** | **978** |

  **AND IT SATURATES.** Twenty times lower buys 370 -> 638; the last five-fold gains 13 %. The
  curve flattens near **~650 mm/a, two thirds of NASA**, so the VALUE is spent -- the same shape
  as the `eps_dry`/`co2_band` sweep. **Two things are fixed at every point and are where the
  missing third lives: the condensate is 78.94 g/m2 in all six runs (the rain never draws the
  cloud down) and `P_conv` is 0.000 in all six against the control's 12.7.** So the residual is
  the dead convective scheme and the missing high cloud, NOT the microphysics constant. Default
  stays 0.5.

  **THE CONVECTIVE GAP IS THE DOWNDRAFT, AND THE SUSPECT NAMED HERE WAS WRONG** (`ATM_MC_DIAG=1`,
  new, print-only, default off, 2026-08-31). It walks the `P_conv` recurrence exactly and charges
  each level only the evaporation it could take, so `G - e_d - e_p = P_conv(ground)` is an
  IDENTITY — it closes against the model's own field to every printed digit (12.7029 vs 12.7029
  control, 0.2867 vs 0.2867 flip). `nm` = 100, gate 0, 24 threads, cos-lat global means mm/a:

  | iter 100 | control | flip |
  |---|---|---|
  | generation `g_p` | 28.28 | **357.79** |
  | `e_d` applied | 15.25 (53.9 %) | **357.51 (99.92 %)** |
  | `e_p` applied | 0.32 (1.1 %) | **0.0000 (0.00 %)** |
  | **ground** | **12.70** | **0.287** |

  **THE FLIP GENERATES 12.7x MORE CONVECTIVE PRECIPITATION THAN THE CONTROL AND THE DOWNDRAFT
  EVAPORATES 99.92 % OF IT.** `P_conv` = 0.287 is not a generation failure — convection fires and
  makes more rain than the shipped model does. `e_p`, sub-cloud evaporation, removes **0.00 %**.
  The drier sub-cloud air is real (RH 0.719 against 0.933) and it acts through the WRONG TERM:
  `e_d` is not a sub-cloud process, it runs at every level from the LFS down, and in BOTH arms
  the flux leaving cloud base already equals the flux at the ground. **All of the loss is
  in-cloud.**

  **AND THE TWO EVAPORATION TERMS ARE NOT WEIGHTED THE SAME WAY** — code-level, not yet measured.
  `e_p` carries the convective area fraction `sigma_p` (~0.045 here) and `e_d` carries nothing, so
  the same water is ~22x more evaporable by the downdraft term than by the sub-cloud one purely
  from a missing area weight; `e_d` has no dependence on `M_d` either, so a downdraft that does
  not exist still evaporates. The test is one knob.

  **THE MISSING HIGH CLOUD, AS AN ICE WATER PATH** (2026-08-31, measured from `output_mcd_ctl`
  and `output_mcd_flip`'s own written slices at iteration 100 — no re-run, two independent cuts):

  | column condensate g/m2 | control | flip | Earth |
  |---|---|---|---|
  | total, 87E cos-lat wt. | 1242 | **69.1** | ~80-120 |
  | total, 28N all longitudes | 1629 | **91.9** | |
  | **LWP** 87E | 1122 | **68.0** | **~50-80** |
  | **IWP** 87E | 119 | **1.15** | **~20-30** |
  | max cell `q_c` | 0.605 g/kg | 0.063 | grid mean over 1x1 deg |
  | levels carrying cloud | 36 of 41 | 20 of 41 | |

  **THE FLIP'S LIQUID IS RIGHT AND ITS ICE IS TWENTY TIMES TOO LOW.** 68 g/m2 of LWP against an
  observed 50-80, with a grid-mean peak of 0.063 g/kg implying a physical in-cloud value once
  divided by `f` ~ 0.2-0.3; the control's 1122 is 15-22x the observation and its 0.605 g/kg peak
  is an ordinary IN-CLOUD number carried as a grid mean. **IWP is the better handle on the
  missing cirrus than the LW forcing is, because it has an observational target**: 1.15 against
  20-30, with 10927 m at 0.00 g/m2 against the control's 9.10.
  And the control's ice is only ~4-6x too high where its liquid is 15-22x, **which is why its LW
  forcing lands at 25.35 W/m2 against Earth's ~25 — bought, not earned**, since `cwp_cap_col`
  then divides the whole column by ~79.

  **WHY NO CIRRUS FORMS ALOFT: FOUR BARRIERS, EACH HIDING THE NEXT, AND THE FIRST IS A DEAD
  STORE** (2026-08-31, `dd35c76`. New knobs `ATM_ICE_COLD`, `ATM_T_FLOOR`, `ATM_RH_MIN`, all
  default = shipped; `CloudFraction::hCrit` flattened; one unguarded fourth root in MLR. Shipped
  branch byte-identical at 1 thread, 20 of 20 files, after every step.)

  1. **THE ATMOSPHERE COULD NOT BE COLDER THAN -37 C.** `initTemperatureData` computed a soft
     asymptotic bound and **overwrote it on the next line** with `t = max(t_00, t_curr)`, so the
     comment promising the asymptote described code that never ran. **22.2 % of the initial air
     column sat at EXACTLY 236.15 K** -- every level from 2987 m up in some columns, 100 % of the
     lid -- and `t_top_init` snapshots the clamped lid so `BC_Atm` re-imposes it every iteration.
     At iteration 100 the coldest cell in the whole atmosphere was still exactly **-37.0000 C in
     every arm**. `t_00` is a PHASE-TRANSITION constant; the physical bound is a tropopause
     temperature (`ATM_T_FLOOR=216.65`).
  2. **ICE WAS DELETED BELOW -37 C IN FIVE PLACES**, one of which also deleted the water vapour
     (`InitValues_Atm.cpp`, `SaturationAdjustment` x3, `UtilsAtm.h`). Cirrus lives at -40 to
     -70 C. Below `t_00` liquid must FREEZE, not vanish; deposition must continue; vapour is not
     a condensate. `ATM_ICE_COLD` fixes all five.
  3. **`H_crit` TURNED BACK UP ABOVE 550 hPa** -- roots at p = 0 AND p = 1000, so 0.50 at 10.9 km
     and 0.68 at the lid, climbing exactly where RH falls. Flattened above the minimum.
  4. **AND THE BINDING CONSTRAINT WAS THE PRESCRIBED HUMIDITY, WHICH THE TEMPERATURE CANNOT
     REACH.** Manabe-Wetherald is linear in sigma and drives RH to ZERO at p -> 0: **0.16 at
     10.9 km against an observed 0.4-0.7 over ice.** Cooling cannot fix it -- the initial vapour
     is `c = RH*q_sat`, so `ATM_T_FLOOR` took 10.9 km from -34.1 C to -42.4 C and moved RH there
     **26.6 % -> 27.1 %**. `ATM_RH_MIN` floors the profile.

  **EACH OF 1-3 IS A NULL ALONE, BECAUSE EACH IS HIDDEN BY THE ONE BEFORE.** The ice cutoff
  cannot bite when the temperature floor is the SAME CONSTANT; the threshold cannot bite when the
  humidity is prescribed below it. `nm` = 100, gate 0, 24 threads, 87E section:

  | arm | LWP | IWP | >7 km | LW forcing | Precip |
  |---|---|---|---|---|---|
  | flip base | 67.98 | 1.15 | 0.035 | 9.02 | 577 |
  | + `H_crit` flat | 67.41 | 1.08 | 0.024 | 8.94 | 586 |
  | + `ATM_ICE_COLD` | 67.86 | 1.28 | 0.026 | 8.91 | 587 |
  | + `ATM_T_FLOOR` | 67.70 | 1.29 | 0.026 | 9.02 | 583 |
  | + `ATM_RH_MIN` 0.35 | 68.26 | 1.62 | 0.463 | 11.02 | 571 |
  | **+ `ATM_RH_MIN` 0.45** | **98.92** | **13.54** | **15.701** | **35.92** | **836** |
  | *Earth* | *50-80* | *20-30* | | *~25* | *978* |

  **ALL FOUR TOGETHER GIVE THIS MODEL CIRRUS FOR THE FIRST TIME**: `q_i` at 10.9 km 0.00000 ->
  **0.00280 g/kg**, condensate above 7 km 0.026 -> **15.70 g/m2**, IWP 1.29 -> **13.54** against
  an observed 20-30, precipitation 583 -> **836 mm/a** against NASA's 978. **0.45 IS PAST THE
  MARK AND IS NOT A CALIBRATION**: the LW forcing overshoots (35.9 against ~25) and the LWP
  leaves its band (98.9 against 50-80). Observed upper-tropospheric RH over ice is 0.4-0.7 and
  the answer is inside that range; do not fit it to the forcing.

  **THE `ATM_RH_MIN` SWEEP, AND IT BRACKETS AT 0.40-0.45** (2026-08-31, `nm` = 100, gate 0,
  24 threads, all four fixes on, 87E section; `ATM_RH_CRIT` = 0.30 throughout):

  | `ATM_RH_MIN` | LWP | IWP | >7 km | RH 7.4 km | RH 10.9 km | LW forc | OLR clear | OLR all-sky | Precip |
  |---|---|---|---|---|---|---|---|---|---|
  | 0 (unfloored) | 67.70 | 1.29 | 0.026 | 32.0 % | 27.1 % | 9.02 | 263.71 | 254.70 | 583 |
  | 0.35 | 68.26 | 1.62 | 0.463 | 35.5 % | 37.2 % | 11.02 | 263.39 | 252.37 | 571 |
  | **0.40** | **74.65** | **4.98** | **5.086** | **40.2 %** | **42.6 %** | **20.48** | **263.21** | **242.73** | **579** |
  | 0.45 | 98.92 | 13.54 | 15.701 | 45.6 % | 48.5 % | 35.92 | 263.03 | 227.11 | 836 |
  | **0.55** | **263.62** | **45.30** | **53.965** | 56.5 % | 60.3 % | **64.72** | 263.09 | **198.38** | **3908** |
  | *Earth* | *50-80* | *20-30* | | *40-70 %* | *40-70 %* | *~25* | *~265* | *~240* | *978* |

  **0.40 IS THE BEST SINGLE POINT AND NOTHING IN THIS TREE HAS EVER PUT THIS MANY QUANTITIES ON
  THEIR OBSERVED VALUES AT ONCE**: clear-sky OLR 263.2 against ~265, **all-sky OLR 242.7 against
  ~240**, LWP 74.7 inside the 50-80 band, and RH 40 / 43 % inside the observed 40-70 %. The LW
  forcing is 20.5 against ~25, low but the closest this model has been.

  **AND THE RESPONSE IS VIOLENTLY SUPER-LINEAR, SO THE USABLE RANGE IS NARROW.** `q_c = f^2*D`
  with `f = (RH - H_crit)/(2(1 - H_crit))` is QUADRATIC in the excess over `H_crit`, and the
  precipitation then compounds it through accretion. IWP goes 1.6 -> 5.0 -> 13.5 -> **45.3** and
  precipitation 571 -> 579 -> 836 -> **3908 mm/a, four times NASA**, for RH_MIN steps of 0.05.
  **0.55 is not "a bit more": it is a different climate.** Do not extrapolate off the top of
  this table.

  **THE THREE TARGETS DISAGREE, AND THAT IS THE NEXT FINDING RATHER THAN A TUNING PROBLEM.**
  Radiation wants 0.40; precipitation wants ~0.45-0.47 (836 at 0.45 against NASA's 978); IWP
  wants higher still (13.5 at 0.45 against 20-30). No single value satisfies all three, and the
  reason is structural: **`ATM_RH_CRIT` is ONE `H_crit` doing two jobs** -- setting the liquid
  deck below 550 hPa and the cirrus above it -- so buying ice costs liquid at a fixed exchange
  rate (LWP 74.7 -> 98.9 -> 263.6 as IWP goes 5.0 -> 13.5 -> 45.3). A separate upper-level
  critical humidity is the obvious next knob, and it is NOT written.

  **THE THRESHOLD SPLIT WORKS, IWP AND PRECIPITATION BOTH LAND, AND IT EXPOSES A NEW DEFECT:
  THE CIRRUS COVERS THE WHOLE GLOBE** (`ATM_RH_CRIT_ICE`, new, default 0 = disabled and verified
  byte-identical over 20 files at 1 thread; `608b3be`). Ramps `H_crit` from `critMid()` at
  550 hPa to `critIce()` at 300 hPa, flat above. All four earlier fixes on, `ATM_RH_MIN` = 0.40:

  | `ATM_RH_CRIT_ICE` | LWP | IWP | >7 km | high-cloud cover | LW forc | OLR all-sky | Precip |
  |---|---|---|---|---|---|---|---|
  | none (0.30) | 74.65 | 4.98 | 5.09 | 58.3 % | 20.48 | 242.7 | 579 |
  | 0.20 | 86.84 | 14.92 | 20.79 | 96.8 % | 42.51 | 220.7 | 738 |
  | **0.15** | 95.47 | **21.39** | 31.89 | 97.9 % | 51.17 | 212.1 | **931** |
  | 0.10 | 105.33 | **29.57** | 45.39 | 100.0 % | 58.61 | 204.8 | 1197 |
  | *Earth* | *50-80* | *20-30* | | *~20-25 %* | *~25* | *~240* | *978* |

  **THE DECOUPLING IS REAL.** RH aloft is 40.6 / 43.3 % at EVERY point -- the knob moves the
  threshold and not the humidity, by construction -- and ice is now much cheaper in liquid: the
  `ATM_RH_MIN` route needed LWP 98.9 to reach IWP 13.5, this reaches IWP 14.9 at LWP 86.8 and
  IWP 21.4 at 95.5. **IWP lands in the observed 20-30 band for the first time, and
  precipitation reaches 931 mm/a against NASA's 978 -- 95 %, from 59 % this morning.**

  **AND THE LW FORCING IS STILL 2x TOO LARGE AT THE RIGHT IWP, WHICH IS A NEW AND CLEANLY
  ISOLATED DEFECT.** 51.2 W/m2 at IWP 21.4 against Earth's ~25 at IWP 20-30. It is NOT the
  optical constants -- `k_liq` = 0.12 and `k_ice` = 0.055 m2/g are already phase-differentiated
  after Stephens (1978). **It is the COVER: the cirrus sits in 97.9 % of columns against Earth's
  ~20-25 %.** The model reaches the right total ice by spreading a thin cirrus over the entire
  globe, and a global thin cirrus forces about twice what a patchy realistic one does.

  **AND THE CAUSE IS THE REPAIR ITSELF, NOT SOMETHING IT UNCOVERED.** `ATM_RH_MIN` is a UNIFORM
  floor: every column is lifted to the same RH aloft, so every column sits the same distance
  above `H_crit` and every column makes the same cirrus. **This is the original "cloud in 99 % of
  columns" pathology reproduced one layer up** -- the shipped model got it from a constant-RH
  column, and this gets it from a constant-RH FLOOR. A flat floor cannot produce patchy cloud.
  The real fix is upper-level humidity carried by the circulation rather than prescribed, which
  is not a knob; until then read `ATM_RH_MIN` as a scaffold and do not tune the LW forcing
  against it.

  **A LATITUDE-DEPENDENT FLOOR, CONFINED TO THE UPPER TROPOSPHERE, AND THE MODEL NOW MAKES
  NASA-MATCHING RAIN FROM A PHYSICALLY SIZED CLOUD** (`ATM_RH_MIN_LAT`, `ATM_RH_MIN_PTOP`, both
  new and default off, each verified byte-identical over 20 files at 1 thread). `ATM_RH_MIN` sets
  the TROPICAL value; the subtropics and storm track follow at the observed ratios 1 : 0.45 :
  0.70 (annual-mean 300 hPa RH over ice ~0.60 / 0.27 / 0.42), via Gaussians at the equator and
  55 deg. `ATM_RH_MIN_PTOP=<hPa>` ramps the floor in over the 200 hPa below it.

  | arm | LWP | IWP | LW forc | OLR all | Precip | 0-15 | 15-35 | 35-65 | 65-90 | global |
  |---|---|---|---|---|---|---|---|---|---|---|
  | flat 0.40, ICE 0.15 | 95.5 | 21.4 | 51.2 | 212 | 931 | 100 % | 100 % | 100 % | 79 % | 97.9 % |
  | lat 0.65, unconfined | 208.0 | 22.1 | 30.9 | 232 | 1558 | 100 % | 18.2 % | 45.7 % | 0 % | 46.4 % |
  | lat 0.65, ICE 0.25 | 216.5 | 28.3 | 37.1 | 226 | 1549 | 100 % | 38.5 % | 82.8 % | 18.4 % | 66.7 % |
  | lat 0.65, PTOP 400 | **82.2** | 14.9 | **29.0** | **234** | 779 | 100 % | **13.0 %** | 45.7 % | 0 % | 44.8 % |
  | **lat 0.65, PTOP 500** | 113.5 | **21.2** | **30.4** | **233** | **1061** | 100 % | **18.2 %** | 45.7 % | 0 % | 46.4 % |
  | *Earth* | *50-80* | *20-30* | *~25* | *~240* | *978* | *~40 %* | *~15 %* | *~35 %* | *~30 %* | *~20-25 %* |

  **THREE SEPARATE LESSONS, IN ORDER.**
  (1) **A FLAT FLOOR PAIRED WITH A LOW ICE THRESHOLD CANNOT MAKE PATCHY CLOUD**, and neither can
  the latitude structure on its own: with `ATM_RH_CRIT_ICE` = 0.15 the threshold sits BELOW every
  latitude's floor, so cover stays 93-98 % however the floor is shaped. The two knobs have to be
  set so the SUBTROPICAL floor falls below `H_crit` and the tropical one well above it -- at the
  unsplit 0.30 the subtropics land at 0.31 and go essentially clear.
  (2) **THE FLOOR WAS NOT AN UPPER-TROPOSPHERIC FLOOR.** Manabe-Wetherald gives RH 0.665 at
  900 hPa and 0.589 at 800, so an unconfined 0.65 BINDS FROM 800 hPa UP -- through the liquid
  deck it was never meant to touch. That is the whole of LWP 208 and precipitation 1558.
  Confining it takes LWP 208 -> 82 and precipitation 1558 -> 779 while KEEPING the cover
  structure and the LW forcing.
  (3) **THE COVER STRUCTURE IS NOW RIGHT WHERE THE CIRCULATION HAS STRUCTURE AND WRONG WHERE IT
  DOES NOT.** Subtropics 18.2 % against ~15 % and mid-latitudes 45.7 % against ~35 % are good;
  the tropics saturate at **100 % against ~40 %** and the poles are **0 % against ~30 %**. The
  tropical failure is the same defect one axis over: **a zonally uniform floor cannot make
  longitudinally patchy cirrus**, and real tropical cirrus is convective and patchy in longitude.
  No prescribed profile can fix that.

  **THE `ATM_RH_MIN_PTOP` SWEEP, AND IT SEPARATES THE TWO PHASES CLEANLY** (lat floor 0.65,
  unsplit threshold, `nm` = 100, 24 threads):

  | `PTOP` hPa | LWP | IWP | LW forc | OLR all | Precip | 15-35 deg cover |
  |---|---|---|---|---|---|---|
  | 400 | **82.2** | 14.85 | 28.97 | 234.4 | 779 | **13.0 %** |
  | 450 | 94.7 | 18.48 | 29.83 | 233.5 | 897 | **15.6 %** |
  | **475** | 103.0 | **20.09** | 30.17 | 233.2 | **992** | 18.2 % |
  | 500 | 113.5 | **21.20** | 30.44 | 233.0 | **1061** | 18.2 % |
  | 550 | 135.2 | 22.02 | 30.74 | 232.7 | 1282 | 18.2 % |
  | unconfined | 208.0 | 22.08 | 30.90 | 232.4 | 1558 | 18.2 % |
  | *Earth* | *50-80* | *20-30* | *~25* | *~240* | *978* | *~15 %* |

  **THE ICE SATURATES BY 500 hPa AND THE LIQUID NEVER DOES.** IWP runs 14.9 -> 18.5 -> 21.2 ->
  22.0 -> 22.1 and the LW forcing 29.0 -> 29.8 -> 30.4 -> 30.7 -> 30.9, both flat from 500
  downward; LWP runs 82 -> 95 -> 113 -> 135 -> **208** and precipitation 779 -> 897 -> 1061 ->
  1282 -> **1558**, both still climbing at the end. **So the cirrus lives above ~500 hPa and
  everything the floor adds below that is liquid** -- which is exactly what the knob was written
  to separate, and it is the evidence that it does. Extending the floor below 500 hPa buys no
  ice and no radiation, only rain.

  **`PTOP` = 475 IS THE PRECIPITATION-OPTIMAL POINT AND IT WAS PREDICTED BEFORE IT WAS RUN.**
  Interpolating 897 (450) and 1061 (500) put NASA's 978 at ~475 with LWP ~104 and IWP ~20;
  measured, it gives **Precip 992.2 mm/a (+1.5 % on NASA), IWP 20.09 (inside the observed
  20-30), LWP 103.0**. That is the first time a value in this chain has been predicted from the
  curve and then confirmed, which is the difference between a sweep and a calibration. **The
  residual is the LIQUID**: even at `PTOP` = 400, where the ice is starved to 14.9, LWP is 82
  against a 50-80 band and the LW forcing 29.0 against ~25. No value of this knob removes that,
  and it is consistent with the tropics still being 100 % covered.

  **THE PRECIPITATION FIGURE IN THE PARAGRAPH BELOW IS PRE-FLIP AND STANDS ON 8129 mm/a OF
MANUFACTURED WATER** (measured 2026-09-01; see *The accepted configuration's 992 mm/a* above).
The same configuration under the 2026-09-01 defaults gives **975.8 mm/a with the floor injecting
0.20** — the radiation, cloud and humidity numbers in it are unaffected, only the precipitation
line is.

**THE BEST CONFIGURATION FOUND (2026-08-31), AND WHAT IT IS WORTH:** `ATM_RH_PROFILE=1
  ATM_RH_CRIT=0.30 ATM_CLOUD_FRAC=1 ATM_CWP_CAP=off ATM_CLOUD_RAD_FRAC=1 ATM_QC_CRIT=0.05
  ATM_ICE_COLD=1 ATM_T_FLOOR=216.65 ATM_RH_MIN_LAT=1 ATM_RH_MIN=0.65 ATM_RH_MIN_PTOP=475`
  gives **Precip 992 mm/a against NASA's 978 (+1.5 %)**, IWP 20.1 in the observed 20-30 band, LWP
  103 against 50-80, cloud LW forcing 30.2 against ~25, all-sky OLR 233 against ~240 and
  clear-sky 263 against ~265. **Every remaining bias has the same sign and about the same size**
  -- LWP +29 %, LW forcing +21 %, precipitable water +21 % (30.3 mm against ~25), all-sky OLR
  -3 % -- which is one excess of low liquid cloud and vapour, not four independent errors, and it
  traces to the tropics still being 100 % covered. **The shipped model matches NASA too -- from a condensate 20x too
  large, cloud in 98.9 % of columns, and a `cwp_cap_col` dividing the column by 79.** The
  difference is that this one does it from a physically sized cloud field. **It is still not a
  prediction**: `ATM_RH_MIN_LAT`'s Gaussians know nothing about where THIS model's ITCZ is, so
  the cover agreement is assumed. Read the whole set as a scaffold that makes the microphysics
  testable, not as a tuned climate.

  **AND LOWERING THE FLOOR EXPOSED A NaN IN THE SHIPPED RADIATION.** `MultiLayerRadiation.h`
  inverts `sigma T^4` as `pow(rad/sigma, 0.25)` **unguarded**, while the identical expression 75
  lines earlier carries `max(1.0, ...)`. The Thomas back-substitution returns a NEGATIVE emission
  on a cold, optically thin column, so the model NaNs at initialisation inside
  `apply_co2_perturbation`, before the time loop. Latent on the shipped branch **only because the
  temperature was clamped at -37 C**. Guarded. This is item 3 again: the solver is not
  energy-closed, so nothing prevents a negative emission.

  **THE TEMPERATURE IS ALSO NOT REACHABLE BY THE RADIATION KNOBS AT THIS `nm`.** `ATM_RAD_EQUIL=1
  ATM_SW_INSOL=1361` as a full-model arm moved the cloud LW forcing 8.94 -> 49.5 W/m2 and left
  the profile **unchanged** (-28.0 vs -28.1 C at 9 km): MLR reaches `t` only through `t_eq` on a
  `teq_refresh_stride` = 20 cadence, and 100 iterations is 20 seconds. The INITIAL profile is the
  only lever on this timescale -- which is what makes `ATM_T_FLOOR` the temperature fix and not
  `ATM_RAD_EQUIL`.

  **RAW SUMS OF `e_d` AND `e_p` ARE NOT A BUDGET.** They are a DEMAND that the `max(0, ...)`
  truncates at every level: summed raw they read 45 000 % of generation with a compensating
  negative clamp, which says only that the demand is large. The unmet remainder is reported
  separately (753x generation under the flip, 2747x in the control).

  *The two figures per point are the two PARITIES*: `printDataAtm` runs every iteration while
  `moist_stride` = 2 runs the ice scheme every other one, so precipitation carries a 2dt
  sawtooth -- **x1.22 in the control, x2.5 under the flip**. Quote both or neither, and do not
  fit a constant against one parity.

  **A NOTE ON HOW THIS WAS FOUND, BECAUSE IT IS THE FAMILY'S OWN PATTERN AGAIN.** The arm was
  first written up here on the condensate and the OLR alone and called a success. The
  precipitation was never looked at -- in a model named ATOM_Precipitation -- and the user found
  it by opening `output_cd400_on` and seeing no rain. Sixth occurrence of "the answer was in
  output nobody opened", and the first where the model's PRIMARY diagnostic was the one skipped.
  **Check `Precip mean` against the NASA line before calling any moisture change an improvement.**


- **THE CONDENSATE IS 20-30x TOO LARGE BECAUSE THE INITIAL HUMIDITY IS NEAR-SATURATED AT EVERY
  LEVEL, AND THE CLOUD SCHEME ONLY WORKS BECAUSE OF IT** (2026-08-29). Traced end to end with
  two new print-only instruments, `ATM_CWP_CENSUS` and `ATM_CLOUD_INIT_DIAG`.

  **The column path and where it sits.** 1584 g/m2 cos-lat mean against an observed ~50-100,
  **98.9 % of columns cloudy**, spread over **38 of 41 levels**, peaking at 127 g/m2 in the
  single level at 3316 m. But the PER-CELL values are ordinary — max cloud water 0.72 g/kg,
  ~0.25 g/kg at the profile peak, max vapour 30.2 g/kg. **Nothing is too wet in any one cell;
  cloud simply exists everywhere.**

  **The humidity is the cause.** `initWaterWapour` sets `c = RH_init*q_sat` with `RH_init`
  **constant in height** (0.60 land / 0.75 ocean) and then multiplies by **1.25** — a fudge whose
  own comment says it "gives a nice cloud around 1 km height". So the ocean column stands at
  **RH = 0.9375 from the surface to the lid.** Measured:

  | z | 82 m | 819 | 1368 | 2163 | 3316 | 4988 |
  |---|---|---|---|---|---|---|
  | mean RH | 0.700 | 0.860 | 0.896 | 0.914 | 0.929 | **0.937** |
  | `H_crit` | 0.983 | 0.934 | 0.899 | 0.859 | 0.823 | 0.801 |
  | cells RH > 0.8 | 70 % | 86 % | 90 % | 93 % | 98 % | **100 %** |
  | cells with cloud | 15 % | 50 % | 90 % | 92 % | 92 % | 78 % |

  **RH RISES with height where Earth's falls** (~80 % in the boundary layer to 40-60 % aloft),
  while `H_crit` FALLS. They cross at ~1.4 km and above that essentially every cell condenses.
  `cloud/cloud_max` is 0.47-0.60 at the peak, so the amount is set by the supersaturation and
  not by the ceiling — the field is a faithful response to the humidity.

  **AND THE TWO CONSTANTS ARE A TUNED PAIR — THE CANCELLING-ERROR PATTERN A THIRD TIME.**
  `ATM_RH_PROFILE=1` (default 0, bit-identical) replaces the constant column with
  Manabe-Wetherald `RH(sigma) = RH_s*(sigma - 0.02)/0.98` and drops the 1.25. Alone it takes the
  column path **1584 -> 0.0006 g/m2**: no cell reaches `H_crit` any more. Lowering the threshold
  with it (`ATM_RH_CRIT`, default 0.8 = shipped) barely helps —

  | `ATM_RH_CRIT` | 0.8 | 0.55 | 0.45 | 0.35 |
  |---|---|---|---|---|
  | column path g/m2 | 0.0006 | 0.0006 | 0.0006 | **1.90** |

  — still ~50x BELOW the target. **The reason is structural and it is the real finding.**
  `cloud_max[i]`, the scheme's ceiling, is itself the horizontal MEAN of `max(0, c - 0.74*q_sat)`
  over the level. In a realistic atmosphere the grid-mean supersaturation is essentially never
  positive, so `cloud_max` collapses to its 1e-4 floor and almost no cloud forms at any
  threshold. **The scheme diagnoses cloud from GRID-MEAN supersaturation, and real cloud forms
  from SUB-GRID variability** — parts of a cell saturated while the mean is not. That is exactly
  what a fractional (Sundqvist/Smith) scheme supplies and what this model has never had.

  **So the chain is: no sub-grid cloud scheme -> the only way to get cloud is to saturate the
  whole column -> the 1.25 fudge does that -> the condensate comes out 20-30x too large ->
  `cwp_cap_col` = 20 divides it by 79 three modules downstream.** Every link is now measured.
  **The repair is a cloud FRACTION, not a re-tune**: nothing in the existing scheme can produce
  a realistic condensate field from a realistic humidity field, so `ATM_RH_PROFILE` must not be
  flipped until the fractional scheme exists. Both knobs stay default OFF and both are
  bit-identical unset (verified: 1583.96 either way).


- **THE TWO MISSING INSTRUMENTS ARE IN, AND BOTH ANSWERED ON FIRST LIGHT** (2026-08-29,
  print-only, verified: physics bit-identical, the only log difference is the new lines).

  **(1) `div(rho u)/rho` is printed alongside `div(u)` at every existing call site.** This file
  said the metric question "is NOT DECIDABLE BY MORE A/B RUNS -- the prerequisite is the
  `div(rho u)/rho` print, not another arm", and the print did not exist. Formed term for term as
  the anelastic Poisson source forms it (`div_src += aux_u * dlnrho[i] * exp_rm`), from the base
  state `densities()` builds, so it prints on BOTH branches and measures what the SHIPPED solver
  leaves behind. First reading, after the initial projection: `div(u)` rms **3.013e-03**,
  `div(rho u)/rho` rms **3.570e-03** — **mass/volume = 1.185**. The mass divergence is 18.5 %
  larger than the volume one, which is the quantity `Psi(ground)`'s non-closure is made of.

  **(2) `div(u)` is now reported EITHER SIDE of `project_velocity_in_loop`.** This file recorded
  that "'the projection zeroes `div(u)` and one time step puts all of it back' and 'the
  projection never reduces `div(u)` at all' have the IDENTICAL signature, and they are different
  defects." A pre/post pair inside one call separates them. At `ATM_PROJECT_IN_LOOP=10`:

  | iter | pre | post | change |
  |---|---|---|---|
  | 1 | 3.013e-03 | 2.748e-03 | **-8.8 %** |
  | 2 | 2.149e-03 | 1.856e-03 | -13.6 % |
  | 3 | 1.951e-03 | 1.675e-03 | -14.1 % |

  **BOTH READINGS WERE WRONG.** The projection does not zero `div(u)` — it removes **9-14 %** of
  it — and it does not fail to act either. And the time step does NOT put it back: post of one
  call is 2.748e-03 while pre of the next is 2.149e-03, i.e. LOWER. **That is why
  `ATM_PROJECT_IN_LOOP` is a null on `Psi(ground)`**: one call dents the divergence by a tenth,
  and `Psi(ground)` is a domain-scale integral of what remains.
  *Caveat*: `nm = 3`, 8 threads, inside the initial transient where the field is still settling,
  so the trend across iterations is spin-up and not a projection result. The pre/post PAIR is
  the clean measurement — it is self-contained within a single call — and it is what settles the
  ambiguity. Re-read the trend on a spun-up run before quoting it.


- **`moist_phys_start_iter` IS 0 SINCE 2026-08-29, AND THE GATE NEVER WORKED.** It was installed
  to defer a MoistConvection/ice-scheme runaway; the recorded blow-up is at **iteration 309 —
  nine iterations AFTER the gate released at 300** (Patagonian Andes, 49S 72W). So the gate did
  not prevent the runaway, it deferred it onto a WORSE field: activating stiff microphysics onto
  a spun-up high-CAPE circulation with a Patagonian jet is what ignited it. It was cured at
  source instead — `MCv_max` 0.5 -> 0.05, the `S_*` and `MC_*` cap passes, and `damp_wiggles` on
  `MC_t`/`MC_v`/`MC_w`. **Validated here before the flip**: `nm = 400`, gate 0, 24 threads —
  400 iterations, **no NaN**, straight through 300-320, `max|w|` **19.83 m/s** against the +-100
  clamp (`max|v|` 2.79, `max|u|` 0.028), `residuum_atm` smooth and bounded 23.55 -> 23.63.
  Set `<moist_phys_start_iter>300</...>` to restore the old branch.
  **What it was silently doing**: with `nm` <= 300 the moist block never executes, and that is
  every A/B config in this tree — so `SaturationAdjustment`, the ice scheme and `MoistConvection`
  have not run in ANY recorded measurement here, `damp_wiggles(t)` included (it sits inside the
  same block).

- **BUT THE GATE IS NOT WHY THE CLOUD CONSTANTS WERE FITTED TO AN INITIAL CONDITION — THE RUN
  LENGTH IS** (`ATM_CWP_CENSUS=1`, new, print-only). Two arms at `nm = 100`, 24 threads, the
  cos-lat-weighted column condensate path in g/m2:

  | arm | mean | p05 | p50 | p95 | max | cloudy (>5 g/m2) |
  |---|---|---|---|---|---|---|
  | gate 300 (no moist physics at all) | 1587.5 | 214.6 | 1967.2 | 2258.3 | 2309.2 | **98.94 %** |
  | gate 0 (moist physics from iter 0) | 1583.3 | 213.6 | 1960.7 | 2252.0 | 2299.1 | **98.94 %** |

  **0.27 % apart.** Turning the moist physics on changes the condensate field by almost nothing,
  because 100 iterations is **20 seconds of physical time** — `SaturationAdjustment` is
  instantaneous and the field arrives already adjusted, while every rate-based process carries
  `dt` and cannot act. This REFUTES the prediction that motivated the census (that arm B would
  differ materially and so `cwp_cap_col` was a fit to the wrong field). The field is the same
  either way; it is `initCloudIce`'s, and the run is too short for anything to move it.

- **AND THE CENSUS FOUND SOMETHING LARGER THAN THE MISSING CLOUD FRACTION.** The column
  condensate path is **1587 g/m2 in 99 % of columns**, against an observed overcast LWP of
  ~100 g/m2 and a global mean nearer 30-50. `cwp_cap_col` = 20 divides that by **79**. So the
  cap is not merely standing in for a cloud fraction, as its own note says — it is compensating
  for a condensate field **20-30x too large covering essentially the whole globe**, and a cloud
  fraction diagnosed from that field would be 1.0 everywhere. The p05 is 214 g/m2: even the
  driest 5 % of columns are twice a realistic overcast column. **The defect is upstream of the
  radiation, in `initCloudIce` and the moist physics, and re-tuning any cloud constant before it
  is fixed is fitting to a field that is wrong by more than the constant.**


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
  **THE HEADING BELOW SAYS THIS WAS FLIPPED ON AND IT WAS NOT** (checked 2026-09-02 against
  `MultiLayerRadiation.h:293` and against the `[RUN CONFIG]` banner of every run in
  `output_rsfull100/`: `RAD_TOPO=0*`). The evidence in it stands; only the flip did not happen.

  **THE ARGUMENT FOR FLIPPING IT ON, 2026-08-28, AND THE LONGER RUN IS WHY — THE 4-ITERATION COMPARISON WAS
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
  **THAT PARAGRAPH WAS WRONG AND THIS ONE REPLACES IT** (2026-08-30, measured). It read: "The
  cause is `q_v` = 7.295 at level 0 and exactly 0.000 in every rock level above. `BC_Atm` Pass 3
  copies the mountain-top humidity down (`c.x[0] = c.x[i_mount]`) while the rock is dry, so in
  `tau_i = tau_dry*dp_i/Sum(dp) + tau_wv*vp_i/Sum(vp)` with `vp_i = c_i*dp_i`, level 0 holds
  essentially the whole sub-surface vapour path and collects a large share of the water-vapour
  optical depth." **It does not. It holds about 1.3 %.**

  `vp_i = c_i*dp_i` and level 0 is the THINNEST layer in pressure -- 4.55 hPa on a Tibet column,
  against 30-64 hPa for the real air above the mountain -- so 7.46 g/kg over 4.55 hPa is a small
  path however dry the rock is. Measured on `output_moist0`'s own written slices at iteration 400,
  as level 0's share of its column:

  | slice | land columns | dry | **water vapour** | condensate |
  |---|---|---|---|---|
  | 87E (`zonal_87`, all latitudes) | 77 | 0.54 % | **1.47 %** (max 3.60) | 0.71 % (max 6.6) |
  | 28N (`longal_62`, all longitudes) | 148 | 0.45 % | **1.27 %** (max 2.78) | 0.28 % (max 2.54) |

  **225 land columns, two independent slices, and the largest water-vapour share anywhere is
  3.6 %.** Nothing here can make `tau_layer` = 7.9.

  **AND ON THE BRANCH WHERE THE 7.9 WAS MEASURED, LEVEL 0 IS NOT IN THE COLUMN AT ALL.** With
  `ATM_RAD_TOPO=1`, `i_mount = i_topography > 0` over land and every loop starts there, so the
  Pass-3 copy cannot reach the radiation by construction. The refuted mechanism was proposed for
  the one branch on which it is structurally impossible.

  **THE REAL CAUSE IS `cwp_cap_col`, AND THIS BULLET ALREADY SAID SO TWO PARAGRAPHS EARLIER.** The
  `tau_dry` note above records it: "The driver is `tau_cloud` ... `cwp_cap_col` caps the COLUMN
  condensate path, and with the rock excluded that same 250 g/m2 is shared among fewer, thinner
  air layers -- same column total, larger per-layer `LWP_i`, and `eps_i = 1-exp(-tau_i)`
  saturates." Two contradictory explanations stood in the same bullet for two days and the
  register carried the wrong one forward as its own item 9. **`BC_Atm` Pass 3 needs no change; it
  is doing its documented job of giving surface-flux code the surface state at level 0.**
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
- **`max epsilon` SITS OVER ANTARCTICA AND TIBET, WHERE THERE IS ALMOST NO CLOUD, AND
  `cwp_cap_col` IS WHY: IT INVERTS THE RANKING OF WHICH COLUMNS ARE OPTICALLY THICK**
  (2026-08-30, measured from `output_moist0`'s own written slices — `nm` = 400, moist gate 0,
  24 threads, no re-run needed).

  The full-field `max epsilon` print puts the maximum at **73S 35E, 2986 m** in 38 of 42 prints
  and at 14S 71W (the Andes) in the other two — high terrain both times, and both far from the
  tropical cloud decks. Read off the 87E cross-section (`0Ma_smooth_Atm_zonal_87_400.vtk`, which
  cuts Tibet) at iteration 400:

  | column | eps-max cell | `q_c` | `q_i` | `q_v` | **Epsilon** |
  |---|---|---|---|---|---|
  | Antarctica 71S | i=23, 2.99 km | 0.011 g/kg | 0.022 | 0.12 g/kg | **0.5927** |
  | Tibet 28N | i=28, 4.99 km | 0.070 | 0.005 | 7.5 | **0.5256** |
  | ocean 10S | i=24, 3.32 km | **0.560** | 0.000 | 11.3 | 0.2487 |

  **The ranking is inverted.** Seventeen times less condensate, twice the emissivity. Within the
  Tibet column it is worse still: level 27 carries 0.244 g/kg at `eps` = 0.0195 while level 28
  carries 0.070 g/kg at `eps` = 0.526 — more cloud, 27x less emissivity.

  **THE INVERTER IS `cwp_cap_col`, AND IT IS A PER-COLUMN NORMALISER.**
  `cloud_scale = cwp_cap_col/cwp_raw` rescales EVERY column to the same 20 g/m2, so a column is
  divided by its own wetness:

  | column | `cwp_raw` | `cloud_scale` | levels w/ condensate | fattest level's share | **path the radiation sees** |
  |---|---|---|---|---|---|
  | Antarctica 71S | 30.1 g/m2 | **0.665** | 9 | 36 % | **7.16 g/m2** |
  | Tibet 28N | 528 | 0.0379 | 16 | 40 % | **8.00** |
  | ocean 10S | 2007 | **0.00997** | 26 | 8.8 % | **1.76** |

  The Antarctic column holds **67x less water** and is divided by 1.5 against the tropics' 100,
  so its thickest layer ends up with **4.1x more optical depth**. A dry polar column also
  concentrates what it has into a few levels at the terrain top (36 % in one) where a tropical
  column spreads it over 26 (8.8 %), which compounds the inversion.

  **AND THE RADIATION IS FAITHFUL TO THE SCALED FIELD, WHICH IS THE PROOF.** Over all 4418
  condensate-bearing cells of that cross-section, `corr(Epsilon, RAW layer path)` = **0.689**
  while `corr(Epsilon, SCALED layer path)` = **0.794**. The scheme is not misreading the cloud;
  the cap has already deleted the relationship between how much cloud a column holds and how much
  the radiation sees.

  For Antarctica the arithmetic closes independently: `tau_layer` = 0.898 there, of which the dry
  share is `dp/Sum(dp)` = 26.4/818.6 -> **0.037** and the vapour term is ~0.0001 (`e_surf` =
  0.157 hPa, so `tau_wv` = 0.003 for the WHOLE column). Essentially all of the 0.898 is
  `tau_cloud`, from 0.033 g/kg of condensate.

  **THIS IS A DIFFERENT MECHANISM FROM THE HIMALAYA `tau_layer` = 7.9 RECORDED BELOW.** That one
  is water vapour — `BC_Atm` Pass 3 copying mountain-top humidity into dry rock — and it is on
  the `ATM_RAD_TOPO=1` branch. This one is condensate, on the shipped `ATM_RAD_TOPO=0` branch, and
  the two were being run together. `q_v` at the Antarctic maximum is 0.12 g/kg; there is no
  vapour artefact to blame.

  **`ATM_RAD_COLDIAG=1`** (new, print-only, default off) prints the max-epsilon cell's four
  optical-depth terms with their percentages, its `q_v`/`q_c`/`q_i`/p/T, its cloud fraction, and
  its column `cwp_raw`/`cloud_scale`. Every number above took slice extraction; from now on it is
  one line in the log. **Fifth instrument-shaped defect in this tree** — `Results_Atm.cpp:40` has
  been printing `max epsilon` and its height all along, and a location is not an attribution.

  **CONSEQUENCE FOR THE FLIP.** `ATM_CWP_CAP` must go off in the same step as `ATM_CLOUD_FRAC` —
  already recorded — but the reason is stronger than "it divides the path by 4". It makes the
  radiation's geography wrong: with the cap on, the optically thickest cell on the planet is over
  the East Antarctic plateau.

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
