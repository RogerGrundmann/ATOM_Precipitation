# ATOM

ATOM (Atmospheric and Ocean Model) is a paleo-climate model implementing a finite-difference Navier-Stokes solver on a spherical shell with RK4 time integration. Vertical coordinate stretching in both models is applied. It includes atmosphere and hydrosphere components with optional turbulence closures (k-ε, k-ω, k-ω SST). Topography and bathymetry are greatfully supported by pygplates. Zero-, One-, Two- and Three-Category-Ice-Schemes (COSMO) by switch are available for precipitation computations. A Saturation-Adjustment-Scheme balances water vapour, cloud water and ice. The ocean model detects in the upper levels where Ekman flow is present upwelling and downwelling regions. The temperature and salinity distributions are controlled as well. In preparation is a deep ocean part as challenge counts here the thermohaline conveyor belt.

## Repository layout

```
atmosphere/   C++ source — atmosphere model
hydrosphere/  C++ source — hydrosphere model
cli/          Command-line interface source and config files
lib/          Shared utilities (arrays, geometry, …)
python/       Python bindings (Cython)
data/         Input data (bathymetry, temperature, precipitation, salinity)
tinyxml2/     Embedded TinyXML-2 XML parser
param.py      Parameter definitions (generates .inc and .pyi at compile time)
Makefile      Top-level build file
```

## Prerequisites

| Dependency | Notes |
|---|---|
| g++ ≥ 7 or clang++ ≥ 6 | C++17 required |
| OpenMP | Typically ships with the compiler |
| Cython | Only needed for the Python interface (`pip install cython`) |
| Python ≥ 3.6 | Only needed for the Python interface |

On Debian/Ubuntu:

```bash
apt-get install g++ cython3
```

On Fedora/RHEL:

```bash
dnf install gcc-c++ python3-cython
```

## Build

```bash
git clone https://github.com/atom-model/ATOM.git
cd ATOM
make all
```

This produces two executables: `cli/atm` (atmosphere) and `cli/hyd` (hydrosphere).

To build only one component:

```bash
make atm
make hyd
```

## Running the CLI

### Atmosphere

```bash
./cli/atm cli/config_atm.xml
```

### Hydrosphere

```bash
./cli/hyd cli/config_hyd.xml
```

### Both coupled

Run the atmosphere first for an initial transfer file, then the hydrosphere, or follow whatever coupling sequence your configuration uses. Each executable reads its own XML config and writes output independently.

Output lands in the `output/` directory organised by time slice.

## Configuration

Both models are configured through XML files. The defaults ship in `cli/config_atm.xml` and `cli/config_hyd.xml`. Copy and edit whichever you need — you only have to include parameters you want to change; everything else falls back to the compiled default.

The top-level structure is:

```xml
<atom>
  <common>
    <time_start>0</time_start>
    <time_end>100</time_end>
    <!-- … -->
  </common>
  <atmosphere>
    <turb_model>k_omega_SST</turb_model>
    <!-- … -->
  </atmosphere>
</atom>
```

Key parameters for the atmosphere:

| Parameter | Example value | Effect |
|---|---|---|
| `time_start` / `time_end` | `0` / `100` | Time slice range (Ma) |
| `n` | `200` | Number of time steps per slice |
| `turb_model` | `k_omega_SST` | Turbulence closure (`k_epsilon`, `k_omega`, `k_omega_SST`) |
| `dr` | `0.003` | Radial grid spacing |

The full parameter set and documentation are in `param.py`. The build step auto-generates C++ and Python stubs from this file.

> **Note:** editing `param.py` regenerates `*Params.h.inc`, which changes the C++ class
> layout. Always `make clean` afterwards — an incremental build links object files with
> mixed old/new layouts and corrupts memory.

## Atmosphere ⇄ ocean SST coupling

Beyond the one-way hand-off (the atmosphere writes a surface wind/temperature transfer
file that the hydrosphere reads as its surface boundary condition), the models can be
coupled **both ways**, so the ocean's own sea-surface temperature (SST) feeds back into
the next atmosphere solve instead of the atmosphere always dictating a frozen SST. Three
pieces make this possible; all default to the classic one-way behaviour and are opt-in.

### 1. Ocean surface thermal boundary condition — `sst_relax_alpha` (hydrosphere)

Controls how strongly the ocean surface temperature is relaxed toward the prescribed
atmospheric SST each iteration:

| `sst_relax_alpha` | Behaviour |
|---|---|
| `1.0` (default) | Hard Dirichlet re-pin — the surface is reset to the atmospheric SST every iteration (the long-standing behaviour, bit-identical) |
| `< 1` | Haney-type flux BC — the surface keeps a fraction of its own tendency, so ocean dynamics can build a real SST anomaly instead of having it erased |
| `0` | Unforced (cold-collapse regime — do not use) |

This is a numerical knob, not a physical restoring timescale (see the note in `param.py`).
A value `< 1` is what lets the ocean carry structure that the coupling below can return.

### 2. Reverse channel — hydrosphere SST → atmosphere `t.x[0]`

The hydrosphere writes its surface SST to an iteration-stamped file
`<slice>_Transfer_Hyd_SST_<iter>.vwtp`. The atmosphere can read it and blend it into its
ocean surface temperature at initialisation:

```
t.x[0] ← (1 − sst_coupling_alpha)·t.x[0] + sst_coupling_alpha·SST_ocean
```

| Parameter (atmosphere) | Default | Effect |
|---|---|---|
| `sst_coupling_alpha` | `0.0` (OFF) | Blend fraction of the ocean SST into the atmospheric surface. `0` = no read (bit-identical to the one-way chain) |
| `hyd_sst_iter` | `-1` | Which hydrosphere SST snapshot to read; `-1` = the latest one present |

The blend is ocean-only, finite-checked, SST-clamped to `[-1.8, 40] °C`, and
**ocean-mean-anchored** (the global ocean-mean is preserved, so the coupling moves SST
*structure* without shifting ocean energy). It is applied before the Held-Suarez
relaxation target is snapshotted, so the returned SST also shapes that target.

### 3. Picard outer loop — `python/run_picard_sst.py`

Wraps the two directions in a fixed-point (Picard) iteration: the initial guess is the
frozen zonal parabola; each round re-runs **both** models and blends the ocean SST back
into the atmosphere until it stops changing.

```
round 0:  ATM (no SST yet) → winds ;  HYD → SST₀
round r:  ATM (blends SSTᵣ₋₁) → winds ;  HYD → SSTᵣ ;  report max|SSTᵣ − SSTᵣ₋₁|
```

```bash
python python/run_picard_sst.py <Ma> --rounds 3 --alpha 0.4 --hyd-relax 0.5 --nm 400
```

| Flag | Default | Meaning |
|---|---|---|
| `--rounds` | `4` | Number of Picard rounds (round 0 = the plain one-way chain) |
| `--alpha` | `0.4` | Atmosphere-side blend fraction (`sst_coupling_alpha`) for rounds ≥ 1 |
| `--hyd-relax` | `0.5` | Hydrosphere `sst_relax_alpha`; **must be `< 1`** or the loop is trivial (the surface is re-pinned and the channel returns the same field) |
| `--anderson` | `0` | Anderson-acceleration history depth `m` (`0` = plain Picard). `m ≈ 2–3` mixes past iterates to reach the fixed point in fewer rounds; engages from round 3 on |
| `--anderson-beta` | `1.0` | Anderson mixing/damping (`1.0` = none; `< 1` adds under-relaxation) |
| `--warm-start` | off | Resume each round's **ocean** from the previous round's restart `.bin` so ocean iterations *accumulate* (`R·nm` effective) instead of re-spinning from the IC every round (async Manabe–Bryan style). Winds/SST target stay fresh; the atmosphere still runs from scratch. Needs `nm` a multiple of `--checkpoint` |
| `--atm-warm-start` | off | **Experimental / A-B only.** Also resume the **atmosphere** from the previous round's `atm_restart` `.bin`, so atm iterations accumulate too (round *r* runs from *r·nm* up to *(r+1)·nm*). Off by default on purpose: unlike the ocean, the surface wind the ocean consumes is best fresh at `nm ≈ 400` and *spins down* with accumulated iterations (trades decay, the westerly-everywhere collapse sets in ~1000→1600), so accumulating past ~1000 total iters can **degrade** the forcing. Use only to measure that effect |
| `--nm` | `400` | Iterations per atmosphere/hydrosphere run |
| `--checkpoint` | `100` | VTK/panorama printout stride, and the convergence-monitor sampling cadence |
| `--tol` | `0.05` | Stop when the fixed-point residual `max|gᵣ − uᵣ|` drops below this (K); equals the round-to-round `max|ΔSST|` in plain-Picard mode |

Each round runs in its own sub-directory. Two different convergences matter: the Picard
**residual** `max|gᵣ − uᵣ|` measures atm↔ocean *self-consistency* (the models agreeing),
while the per-round `convergence.csv` (report-only ocean-mean-T and kinetic-energy drift)
measures *physical steadiness* (each model reaching its own equilibrium). The ocean is
slow, so at moderate `nm` the residual can be small while the ocean is still spinning up —
`--warm-start` is the lever that lets ocean iterations accumulate toward equilibrium.

## Correctness fixes ported from ATHAD (2026-08-26)

ATHAD (the Hadean fork of this model) found five defects that are **live here and were latent
on Earth** — three data races and two memory-safety defects — and verifying them here
turned up a sixth, in the ocean's own min/max reporter. They are fixed in this tree as of 2026-08-26.

**EVERY NUMBER RECORDED IN THIS README OR IN A RUN LOG BEFORE THAT DATE WAS PRODUCED BY A
RACY BINARY**, and a re-run will not reproduce it. The races themselves moved the last digits;
**the cure moves more than that**, because red-black is a different sweep order from
lexicographic Gauss-Seidel and converges to the same solution by a different path. Measured at
4 iterations, 4 threads, before against after: `residuum_atm` 21.48903401 -> 21.48908287
(2.3e-6), but **max u-component 1.310076 -> 1.321737 m/s, 0.89 %**, at the same cell. Nothing
here is retracted -- the old runs were not wrong, they were not repeatable -- but a figure
quoted to better than ~1 % from before this date should be re-measured, not cited.

| # | Site | Defect |
|---|---|---|
| 1 | `atmosphere/PressureSolverAtm.h`, `hydrosphere/PressureSolverHyd.h` | The Poisson relaxation wrote `p_dyn` **in place** under `collapse(2) schedule(dynamic, 4)` over (i,j) — the two indices its own stencil reads across. Cured with a red-black colouring of (i+j+k); `schedule(static)`. The family's third encounter with this defect (ATURAN `ffd0e0e` cured it in the shared `PressureSolver.h` by serialising, which is affordable there and is not here) |
| 2 | `atmosphere/UtilsAtm.h`, `hydrosphere/UtilsHyd.h` | `residuum_old` — a shared model member — was written from inside every thread's loop. The reported "old" residuum was whichever thread wrote last, and it drives the *declining / too high* line a human reads to judge convergence. Now captured at entry and written once after the reduction. The reported error **location** got a deterministic (i,j,k) tie-break at the same time |
| 3 | `lib/Utils.cpp` | `m_node_weights` was built lazily inside `GetMean_2D/3D`, which `printDataAtm()` calls from ~9 concurrent `omp sections`. Latent here **by accident** — an earlier single-threaded `GetMean_2D(temperature_NASA)` in `initTemperatureData` happens to populate the table first. ATHAD has no NASA field, so `printDataAtm` became the first caller and the race went live (heap corruption). Now a function-local static, which C++11 guarantees is initialised exactly once |
| 4 | `atmosphere/MoistConvection.h` | `findCloudBaseLFS` indexed `m.t.x[-1][j][k]` when neither pressure scan fired, leaving the index at its `-1` sentinel. Latent here because `p_surf ~1013 hPa` always trips the 1000 hPa threshold by level 1; live in ATHAD, whose column sits entirely above the absolute-hPa thresholds |
| 5 | `atmosphere/InitValues_Atm.cpp`, `hydrosphere/FileIO_Hyd.cpp` | `get_temperatures_from_curve` dereferenced `m.begin()` and decremented `m.end()` **before** the `m.size() < 2` guard — undefined behaviour on an empty map. The guard written to catch a too-small map could not run until after the code it was guarding. The order is now the other way round |

| 6 | `hydrosphere/MinMax_Hyd.cpp` | `searchMinMax_2D/3D` reduce over threads with `#pragma omp critical` and a plain `>`, so the first thread into the section won an equal extremum and **the reported location of a tied max/min depended on thread arrival order**. Found by the smoke test written to verify #1-#5, not by inspection. Ties are what a cap or a floor produces, so the fields most likely to tie are the ones worth watching: `min water density = 997.000000` (the floor), `min salinity = 0.000000`, `max temperature = 0.000000`. Now broken toward the lexicographically smallest (i,j,k), in the per-thread scan and in the merge. The atmosphere's `searchMinMax_*` is a serial loop and never had this |

And one build fix that is the reason the others were worth doing:

- **`Makefile` now passes `-MMD -MP` and `-include`s the generated `.d` files.** Without them
  `make` saw only `.cpp` timestamps, and nearly all the physics in this tree lives in headers —
  so a header edit produced a link of stale objects that looked like a successful build.

  **AND THAT REPAIR HAD A HOLE, WHICH BIT ON THE SAME DAY** (2026-08-26). Header tracking only
  covers objects that have already been compiled *with* `-MMD`; an object built before it has
  no `.d`, so `make` still sees only its `.cpp` and never rebuilds it. `cli/atm.o` was in
  exactly that state — dated 2026-07-28 through every build of the port — and it is the one
  object that matters most, because it holds `cAtmosphereModel model;` as a **stack local**.
  Adding six diagnostic `Array` members to the class made `main` reserve the OLD `sizeof`
  while `libatom.a` constructed the NEW one: the constructor ran off the end of `main`'s frame
  and tripped the stack canary at exit — `*** stack smashing detected ***`, after a complete
  and apparently successful run, with no compiler warning anywhere. Diagnosed by
  disassembling `main` (object at `rbp-0x1600`, canary at `rbp-0x18`, object 0x1678 long).
  AddressSanitizer found nothing, because an object overrunning its own stack slot by ODR
  mismatch is not an out-of-bounds access to anything ASAN tracks.

  The durable fix is one line: **every object now depends on the `Makefile` itself**, so a
  build-flag change forces the full rebuild that generates the missing `.d` files. Worth
  porting to ATHAD, ATHAD_COND and ATHAD_PERID, which have `-MMD` but not this.

**Verification, measured on this tree rather than argued from ATHAD's.** A worktree of the
pre-fix `HEAD` was built beside the fixed tree and each ran the same 4-iteration case twice at
4 threads:

| | pre-fix, run A | pre-fix, run B | fixed, run A | fixed, run B |
|---|---|---|---|---|
| `residuum_atm` (1st call) | 21.48903401 | **21.48903394** | 21.48908287 | 21.48908287 |
| `residuum_old` (1st call) | 19.33304066 | **19.33304053** | **0.00001000** | 0.00001000 |
| max u-component | 1.310076 m/s | **1.310122 m/s** | 1.321737 m/s | 1.321737 m/s |

Read the `residuum_old` row first: on the FIRST call there is no previous residuum, so the
correct value is the `1.0e-5` initialiser. The racy binary printed 19.3 — a value from inside
the call it was supposed to precede, and a different one each run. The `max u-component` row is
the one that matters for results: **a printed physical quantity differed run to run on the same
binary at the same thread count**, by 4.6e-5 m/s.

After the fix the two runs are **bit-identical** (only wall-clock timings differ), and a 2-thread
run gives the identical residuum series as well. What remains is floating-point reduction order,
inherited and not addressed here: a tropical precipitation probe differs 4e-6 relative between 2
and 4 threads (`P_rain[i_check]` 15.4734034 against 15.4734603). ATHAD documents the same residue
and traces its amplification to the convective triggers.

## Instruments ported from ATHAD (2026-08-26)

Three diagnostics, all **print-only or plot-only** — verified: with them in, the residuum
series and every printed extremum are identical to the binary without them. What they measure
is not.

### `checkRadialMetric()` — is `exp_rm` the Jacobian of the radial stretch?

It is documented as one in `TurbulenceAtm.h` and written as one in `PressureSolverAtm.h`
(`dp/dr_physical = exp_rm * dp/d(rad.z)`), and it is not: `exp_rm = 1/(rm+1)` is the Jacobian
of a **quadratic** stretch while `init_layer_heights` builds an **exponential** one. The test
is unit-free — `[inv_2dr * exp_rm] * (z[i+1] - z[i-1])` must be constant in `i` — so it cannot
be argued away as a units convention. Printed at every startup:

    AGCM: radial metric check - core length unit runs 806 m at the surface to 18717 m at
    the top, spread 23.21x (1.00 = exp_rm is the Jacobian).

**23.21x, and that is the worst in the family** — ATHAD is 11.8x at `zeta = 3.0`, ATHAD_COND
~12x — because this tree's `zeta = 3.715` is the largest stretch anywhere in it. Every radial
derivative in the core is mis-scaled by a factor varying 23x across the column, low at the
surface. It is **not** a function of `im`. `ATM_METRIC_CHECK=1` adds the per-level table; the
whole diagnostic goes quiet by itself once the spread drops below 1.05, so it costs one startup
line and cannot change a result. See ATHAD README item 39 for the measurement and for why
`zeta` rather than `im` is the lever.

### The meridional streamfunction is density-weighted

`MinMax_Atm.cpp` applied one constant density — `r_air`, the SURFACE value — at every level of
a mass-flux integral, i.e. it labelled a volume flux kg/s. The density is inside the integral
and inside the zonal mean now, which also keeps the `<rho'v'>` correlation term a warm rising
branch carries. Measured at 4 iterations:

| | before | after |
|---|---|---|
| `Psi_max` | 851.68 (1e9 kg/s) | **373.84** |
| at | lat 45.00, z = 1729 m | lat 45.00, z = 1729 m |

**2.28x overweight at the cell core**, with the location unmoved — so here it distorted cell
STRENGTH rather than hiding a cell, which is what it did in ATHAD's 250 bar column where rho
spans four orders of magnitude instead of six-fold. This tree is the one that uses this
diagnostic to judge exactly that question (`24ff23a` "revive Hadley/Ferrel cells", `1e59daa`
jet spin-down), so **every Psi number recorded here before 2026-08-26 is the old,
volume-flux-like quantity and is not comparable with what the run log prints now.**

### The zonal VTK slice is drawn on a height axis

`paraview_vtk_zonal` wrote the vertical coordinate as **level index**, which on this
exponentially stretched grid stretches the bottom and squashes the top by ~34x, varying with
altitude — so no ParaView aspect setting could undo it, and contours, glyph angles and
streamline curvature all inherited it. It is true height now, mapped onto the same plot span,
so the vertical exaggeration is one constant.

Two unit defects went with it. The glyph vectors were scaled by `1/u_0` while every scalar
written beside them used `u_0` — a factor of `u_0^2` between a component and its own vector, in
the same file (direction unaffected, magnitude not). And raw m/s glyphs on a strongly
anisotropic flow all lie flat along the latitude axis, so a new `uv_plot` vector carries the
field in plot-units per day, in which a closed cell is drawn closed and glyphs agree with
streamlines by construction. `u-v-Cell` is kept, now correctly dimensional.

**The longitudinal writer is still on index space** in this tree and in ATHAD. ATHAD fixed the
zonal one only, and this port does not go beyond what was measured there.

## The radial momentum budget (`ubud_*`), ported 2026-08-26

`vbud_*` and `wbud_*` decompose the meridional and zonal momentum equations; the RADIAL one
had no instrument in either tree until ATHAD's README item 42, where a spurious radial
acceleration of rms 293 went unseen for exactly that reason. Six arrays, captured on
checkpoint iterations only, summing to `rhs_u`: `ubud_advv + ubud_advh` is `-transport_u`
exactly, so the split can be checked rather than trusted. Print-only.

Measured here at 4 iterations, 8 threads (domain maxima of |term|, nondimensional):

    ubud_pgf   = 0.000e+00
    ubud_cor   = 0.000e+00
    ubud_advv  = 1.144e-03   (0.216 of the largest)
    ubud_advh  = 3.058e-03   (0.577)
    ubud_diff  = 5.304e-03   (1.000)
    ubud_buoy  = 1.464e-07   (2.8e-05)

**`ubud_cor` being identically zero is the check the split exists for**: the non-traditional
Coriolis is genuinely off, not merely documented off. ATHAD gets the same.

**`ubud_buoy` answers the question `ATM_BUOY_TREF` was waiting on.** The buoyancy body force is
**2.8e-5 of the dominant term**, and even scaled to full strength (`buoyancy_ramp` is 0.013 at
iteration 4) it would be ~2e-3 of the diffusion — 0.2 %. A **5 % correction to that term is
worth ~0.01 % of the radial momentum balance**, so the 300-iteration A/B item 11 was waiting
for is not worth its wall clock. That is what the budget was ported to decide, and it decided
it in one 40-second run instead of fifty minutes. ATHAD reached the same conclusion from the
other side: `ubud_buoy` ~1e-6 against a pressure gradient of 1.15.

**One caveat on `ubud_pgf = 0`, so it is not over-read.** The capture happens inside RK4, which
runs before the pressure solver in each iteration, and at iterations 2 and 4 `p_dyn` is still
identically zero at that moment (it reaches 0.117 hPa by the end of iteration 4). This is a
spin-up snapshot, not a claim that the radial pressure gradient is absent — in ATHAD it is the
dominant term. Only the buoyancy/diffusion ratio should be read from the table above.

## Tier C knobs ported from ATHAD (2026-08-26)

These change physics, so all of them are **default-off / default-1 and bit-identical when
unset**, and each needs its own A/B *here* rather than inheriting ATHAD's default. The first
one already shows why.

### `ATM_PRESS_SWEEPS`, `ATM_PROJ_SWEEPS` — and ATHAD's lever is INERT in this tree

The initial projection calls the solver 200 times; `ATM_PROJ_SWEEPS` sets the relaxation sweeps
*per pass*, and `ATM_PRESS_SWEEPS` does the same for the time loop. They are deliberately
separate: one knob would have made "10 sweeps in the time loop" also mean 2000 relaxations at
startup, so the arms of any comparison would differ in their initial state as well as in the
quantity under test. ATHAD lost an attribution exactly that way.

The meridional streamfunction must vanish at the ground — `u` is zero at the surface and Psi is
zero at the lid, so the column-integrated meridional mass flux is forced to vanish. It does not,
here or in ATHAD. Measured at 4 iterations, 8 threads (RMS over latitude, not max — the max sits
inside one cell and hides changes elsewhere):

| `ATM_PROJ_SWEEPS` | RMS Psi(ground) | vs 1 | projection cost |
|---|---|---|---|
| 1 (default) | 1.5841e11 kg/s | — | 2.6 s |
| 10 | 1.5834e11 | **-0.04 %** | 21.0 s |
| 100 | 1.5830e11 | **-0.07 %** | 223 s |

Each arm wrote to its own output directory; the 100-arm was re-run alone after an earlier
attempt shared one with a concurrent run, and reproduced 1.5830e11 exactly.

**In ATHAD the same knob cuts it 52.5 %.** The knob is connected — the projection cost rises 8x
and 96x — so this is a real null, and the default stays at **1** here. What it means is that
this tree's Psi(ground), **42 % of its own interior maximum**, is not a projection-convergence
problem at all: it is the structural non-closure ATHAD's item 72 attributes to the discrete
divergence and gradient not being adjoint on a collocated stencil, with Rhie-Chow face
reconstruction named as the un-done repair. In ATHAD roughly half the non-closure was
convergence and the rest structural; here it looks like all of it.

### `ATM_GRID_PRESSURE` — ported with two REFITTED constants, and it has nothing to give

Levels placed uniformly in `ln p` on a reference hydrostatic column instead of exponentially in
height. **Default off, unset bit-identical.** Two constants had to be refitted for this tree,
and the first would have been a silent disaster if copied:

| | ATHAD ships | this tree |
|---|---|---|
| `ATM_GRID_PTOP` | 1e-6 | **0.08538** |
| `ATM_GRID_BETA` | 4.33 | **3.988** |

`p_top = 1e-6 p_0` on Earth's column sits at **81.9 km**. Importing ATHAD's default would not
have regridded this model, it would have **extended the shell 5.1x** and made it a different
model while looking like a grid option. 0.08538 is the pressure at this tree's own legacy lid
(16023 m, 86.5 hPa). `beta = 3.988` is then fitted so the bottom layer matches the legacy
38.9 m; at `beta = zeta = 3.715` it comes out 1.23x coarser.

`buildReferenceColumn` also had to be rewritten rather than copied. ATHAD integrates a **dry
adiabat** because its column is one by construction; on Earth that law puts the 16 km lid at
117 K, and the ladder would be built for an atmosphere this model does not have. It integrates
this tree's own initialisation law instead — the 6.5 K/km reference lapse of
`InitValues_Atm.cpp`, isothermal above the mean tropopause.

**THE OUTCOME WAS PREDICTED ANALYTICALLY BEFORE ANY CODE WAS WRITTEN, AND THE PREDICTION HELD.**
For an isothermal column a `ln p` ladder IS a height ladder; the two differ only through the
temperature contrast up the column, and this tree has the mildest in the family. The ladder has
**Lambda = 2.46 e-foldings** to redistribute across the shell against ATHAD's 13.8, because
**shell/H = 2.01** here (ATHAD 5.06, ATHAD_COND 7.74, ATHAD_PERID 15.54 — and PERID is the one
tree where this branch paid).

Measured at 40 iterations, 8 threads, against the legacy grid:

| | legacy | pressure grid | change |
|---|---|---|---|
| closure ratio | 0.4514 | 0.4570 | **+1.2 % WORSE** |
| max abs Psi(ground) | 3.9454e11 | 3.9968e11 | +1.3 % |
| max abs Psi above 2 km | 3.6363e11 | 3.5678e11 | **-1.9 %** (weaker circulation) |
| `residuum_atm` | 22.3615 | 22.4690 | +0.48 % |
| mean T / KE | 254.379 K / 36.160 | 254.338 / 36.707 | -0.04 K / +1.5 % |
| max radial wind | 0.061202 m/s | 0.060951 | -0.4 % |
| `checkRadialMetric` spread | 23.21x | **21.87x** | -5.8 % |
| bottom layer / lid | 38.9 m / 16023 m | matched / 16023 m | by construction |

**About one percent, and the sign is against it**: the streamfunction closes slightly worse and
the interior circulation is slightly weaker. **ATHAD's free side-benefit does not appear either**
— there the ln-p grid cut the metric spread 11.77x -> 2.19x, a 5.4x improvement that item 82
called a third route to item 39's open decision; here 23.21x -> 21.87x, which is nothing.

So the branch is available, correct, and **not recommended for this tree**. That is a result
about geometry rather than about tuning: a shell two scale heights deep has almost no mass
redistribution available to it, and no choice of `beta` changes that — `beta` moves where the
levels sit within the ladder, and the ladder itself is nearly the legacy one.

### `ATM_METRIC_EXACT` — and ATHAD's result does NOT reproduce here

`exp_rm = 1/(rm+1)` is documented as the Jacobian of the radial stretch and is the Jacobian of a
QUADRATIC one, while `init_layer_heights` builds an EXPONENTIAL one — the defect
`checkRadialMetric` reports as a **23.21x** spread above. `ATM_METRIC_EXACT=1` replaces it with
the true `metricShellLength()/J`, `J = zeta*L_atm*exp(zeta*(r-r0))`, at all twelve sites that
computed it. **Default off; unset is bit-identical** (verified: twelve rerouted call sites,
eleven modified diffusion terms and a new Poisson stencil term, and the only difference against
the pre-refactor binary is a timing line).

**A second defect comes with the first, exactly as in ATHAD.** The radial Laplacian on a
stretched grid is `exp_2_rm*(f'' - (J'/J)*f')` and the core computes `f'' * exp_2_rm` and stops,
so the curvature term is missing under BOTH metrics — small under the legacy one
(`J'/J = 1/(rm+1) <= 0.5`), the same order as the retained term under the true one
(`J'/J = zeta = 3.715`). It is added at the eleven diffusion terms and as a first-derivative
off-diagonal in the Poisson operator. `metricCurv()` returns 0 on the legacy branch, so the
missing legacy term is recorded rather than silently fixed.

With the knob on, `checkRadialMetric` reports the core length unit as **16046 m at every level,
spread 1.00x** — the unit-free test passing by construction.

**Measured at 40 iterations, 8 threads, each arm in its own output directory:**

| | legacy | exact | change |
|---|---|---|---|
| RMS Psi(ground) | 1.6414e11 kg/s | 1.6404e11 | **-0.06 %** |
| max abs Psi above 2 km | 3.6363e11 | 3.6360e11 | -0.008 % |
| closure ratio | 0.4514 | 0.4512 | -0.04 % |
| `residuum_atm` | 22.36153 | 22.34504 | -0.07 % |
| mean T / KE | 254.379 K / 36.160 | 254.378 / 36.159 | ~0 |
| **max radial wind** | **0.061202 m/s** | **0.019108 m/s** | **-69 %** |

**ATHAD's item 80 finding does not reproduce in this tree.** There, making the metric exact left
the OLR untouched but drove `div(rho u)/rho` rms from 2.739e-02 to 7.722e-02 and threw `Psi_max`
from 36 km to the ground — 2.8x worse, and the reason its default is off. Here the same repair
is a **null on every integrated quantity, and marginally in the RIGHT direction**: the residuum
and the closure ratio both fall slightly. Two trees, measured the same way, opposite outcomes —
so *"the correct Jacobian makes the projection worse"* is an ATHAD property and not a family
law. This tree has the WORSE metric (23.2x against 11.8x) and pays nothing to fix it.

The one thing that does move is the **radial-wind extremum, -69 %**. Read it as a spike being
suppressed rather than as a change in the circulation: it is a max over the domain, it sits at a
different cell in the two arms, and everything integrated is unmoved at the 1e-3 level. The
radial component is exactly what a radial metric rescales, so a localised vertical spike
responding while the mass circulation does not is the expected shape.

**The default stays OFF, and the reason is run length, not the result.** Forty iterations is a
spin-up: this tree's own history is full of differences that appear later, and ATHAD's item 54
watched a 4.4 % gap grow to 9.1 % between iterations 20 and 40. What can be said is that the
principled metric is CHEAP here, which makes this the natural tree to test it in properly —
the opposite of the situation in ATHAD.

**Note also what is missing to judge it well.** ATHAD reads this question through
`div(rho u)/rho`, printed every iteration; this tree has no such diagnostic, so the closure
above is the Psi(ground) RMS computed from the CSV. Porting that print is the obvious
follow-on if the question is pursued.

### `ATM_BUOY_TREF`, `ATM_BUOY_CONSISTENT` — implemented, UNMEASURED in this tree

The Boussinesq buoyancy divides its temperature anomaly by `t_0` = 273.15 K, the
**non-dimensionalisation constant**, where the physical reference temperature belongs.
`ATM_BUOY_TREF` divides by `t_ref_level[i]` instead. Here that is a **~5 %** correction
(`T_ref/t_0` = 1.055 at 288 K); in ATHAD it runs 5.49x at the surface and 0.96x at the top, so
there it is a height-dependent distortion rather than a rescaling.

`ATM_BUOY_CONSISTENT` additionally removes the second of item 34's extra `*dt` factors,
`g*dt/u_0 -> g*L_atm/u_0^2`. **That is a factor of 5.0e5 on a body force in this tree.** The
largest coefficient ever run on this term in the family was an intermediate 336, which drove a
polar vertical runaway; expect instability, and if it appears, that is a Boussinesq result
rather than a fault in that line — do not tune the coefficient back down, which is how the 336
arrived.

**Why no measurement yet, stated so it is not mistaken for a null.** A 4-iteration A/B shows
`ATM_BUOY_TREF` changing nothing, and that test is under-powered rather than informative:
`buoyancy_ramp` is **0.013** at iteration 4 (it reaches 1 at `buoyancy_ramp_iters` ~ 300), so
the whole term is suppressed ~100x during the window measured. A real test needs ~300+
iterations per arm, and the *better* instrument is a radial momentum budget — ATHAD's `ubud_*`,
which exists in neither tree until it is ported — because it answers the prior question of how
large this term is at all. ATHAD measured its own as ~1e-6 against a pressure gradient of 1.15,
i.e. effectively absent, on the same shipped coefficient this tree uses.

## Atmosphere diagnostics and A/B knobs

A few atmosphere-model behaviours can be probed at run time without recompiling. These are
research/debugging aids aimed at the circulation spin-down problem (the trades and the
extratropical jet weakening over a run). All are **off / bit-identical by default** — set the
environment variable before launching the atmosphere to change behaviour.

### Environment A/B knobs

| Env variable | Default | Effect |
|---|---|---|
| `ATM_METRIC_RADIUS` | `r_Earth` (on) | Planetary radius, in km, for the HORIZONTAL metric. Since 2026-07-28 the default is the configured `r_Earth`; `0` restores the old `rad.z ~ 1.5` metric for A/B work. It replaced `ATM_CORIOLIS_SCALE`, which was a multiplier compensating for the same length error: that error factorised as 397.5 (the metric, fixed by this knob) x 40 (`force_nd` dividing by `L_atm` = 400 m instead of the ~16023 m one `rad.z` unit represents, now fixed in RHS_Atm_Turb). With both corrected there is nothing left to sweep |
| `ATM_RADIAL_SHAPIRO_STRENGTH` | `1.0` | Scales the strength of the per-iteration radial (vertical) Shapiro filter applied to `u, v, w`. The column-integrated momentum budget identifies these passes as the dominant net sink of extratropical-jet momentum. Values `< 1` ease the filter to test whether that slows the jet spin-down; `0` disables it entirely (**risks** the radial 2Δ checkerboard / near-surface CFL blow-up the filter guards against) |

Example — run the atmosphere with a stronger Coriolis and a gentler radial filter:

```bash
ATM_RADIAL_SHAPIRO_STRENGTH=0.5 ./cli/atm cli/config_atm.xml
```

Each knob is read once (first use) via `getenv`, so it applies to the whole run.

### Momentum-budget diagnostics

Every `checkpoint` iterations the atmosphere writes zonal-mean momentum budgets that
attribute each per-iteration velocity change to a specific term, so a spin-down can be
pinned to a source (e.g. weakening Coriolis) or a sink (a specific filter / diffusion):

- `v_momentum_budget_<iter>.csv` — **meridional** wind (Hadley/Ferrel branches)
- `w_momentum_budget_<iter>.csv` — **zonal** wind (the trades / Walker component)

Each row is one `(latitude, height)` cell. The columns split the net Δ into the RK4 dynamics
(`dw_dyn`) versus each post-RK4 filter (`dw_polar`, `dw_orog`, `dw_radial`), and further break
the RK4 net into its physical terms (`pgf`, `coriolis`, `adv_vert`, `adv_horiz`, `diffusion`,
`drag_mc`) in m/s per iteration. A one-line tropical trade-layer summary is also echoed to
stdout each time.

### ParaView field dumps

The atmosphere's ParaView writers (`atmosphere/Paraview_Atm.cpp`) expose many optional fields
as commented-out `dump_*` lines — force components (`BuoyancyForce`, `CoriolisForce`,
`CentrifugalForce`, `PresGradForce`), `CO2-Concentration`, latent/sensible/radiative heat, and
more. Uncomment the ones you need for a given investigation and rebuild. The pole-singular
`paraview_sphere_vts` output can be skipped in favour of the `panorama` `.vts` (whose
derivatives are well-behaved away from the poles).

## Output

Each run writes output files to `output/`:

- **VTK files** — radial, zonal, and longitudinal slices readable by [ParaView](https://www.paraview.org/)
- **XYZ grid files** — tab-separated gridded data for post-processing

## Python interface

Build the Cython bindings (this also rebuilds `libatom.a`):

```bash
make python
```

This produces `python/pyatom.<platform>.so` in place. Import it from the `python/`
directory (or add that directory to `PYTHONPATH`):

```python
from pyatom import Atmosphere, Hydrosphere

atm = Atmosphere()                       # no-argument constructor
atm.load_config("config_atm.xml")        # optional: load an XML config
atm.run()                                # run every time slice in the config
```

The models can also be configured directly through properties and driven one paleo
time slice at a time — this is what the run scripts do:

```python
atm = Atmosphere()
atm.nm = 400                             # iterations for this slice
atm.checkpoint = 100                     # VTK/MinMax printout stride
atm.output_path = b"output_100Ma/"       # note: paths are bytes in Python 3
atm.run_time_slice(100)                  # run the Ma = 100 slice
```

Every parameter in `param.py` is exposed as a read/write property of the same name.
`run()` iterates all slices defined by the config; `run_time_slice(Ma)` runs a single
slice. `Hydrosphere` has the same interface (plus an `input_path` for the atmosphere
transfer files it reads).

Example drivers in `python/`:

| Script | What it does |
|---|---|
| `run_atm_ma.py` / `run_hyd_ma.py` | Run a single atmosphere / hydrosphere paleo slice |
| `run_chain.sh` / `run_paleo_chain.sh` | Atmosphere → hydrosphere one-way chain |
| `run_picard_sst.py` | Two-way SST Picard loop (see the coupling section above) |

## Input data

The `data/` directory ships with:

- Paleotopography/bathymetry grids 0–140 Ma (Smith et al. 1994; Golonka et al. 1997)
- Present-day surface temperature (NASA)
- Present-day precipitation (NASA)
- Present-day salinity (NASA)

## Authors

Roger Grundmann, Michael Chin
