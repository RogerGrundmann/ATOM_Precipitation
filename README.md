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

## Three more instruments, and what they said (2026-08-26)

All print-only or field-only. `brunt_N2`, `tau_above` and `tau_layer` are dumped by all four
ParaView writers (panorama, zonal, longitudinal, radial) and printed by `print_min_max_atm`;
`div(u)` prints at two call sites in the pressure solver.

### `brunt_N2` — the first ABSOLUTE check this tree has ever had

`N^2 = (g/theta) d(theta)/dz`, potential temperature from `T*(p_0/p)^kappa` with
`kappa = R_Air/cp_l`. **The vertical derivative uses `get_layer_height()`, not the core's
`exp_rm`** — deliberately, so the field is physical and can be judged against reality rather
than inheriting the 23.2x metric error `checkRadialMetric` reports.

Every other test in this README compares the model against itself. This one has an external
answer: Earth's troposphere is near 1e-4 s^-2, its stratosphere near 4e-4.

    level   height    median N^2      statically unstable
        0        0    -3.20e-05           54.1 %      superadiabatic surface layer
        6      298     1.544e-04           0.0 %
       12      819     1.542e-04           0.0 %
       24     3316     1.537e-04           7.2 %
       30     6088     1.695e-04           0.0 %
       33     8173     2.867e-04           0.0 %
       36    10927     4.383e-04           0.0 %
       39    14567     5.161e-04           0.0 %

    whole slice: median 1.547e-04,  3.5 % statically unstable

**It passes.** The troposphere sits at 1.54e-4, flat from 300 m to 6 km; the tropopause appears
without being told where it is, N^2 climbing through 2.9e-4 at 8 km to **4.4e-4 at 10.9 km**,
which is the textbook stratospheric value at the textbook height; and the 3.5 % of statically
unstable air is concentrated at the surface, where a superadiabatic layer belongs.

**Read the printed min/max as boundary-layer values**, not as the field: they are ±0.03-0.047
at 0 m every time. The interior is the signal.

### `tau_above`, `tau_layer` — exposure, not new physics

`MultiLayerRadiation` already computed the per-layer optical depth and threw it away after
converting to `epsilon`. It cannot be recovered afterwards, because `epsilon` saturates at
`tau ~ 37`, which is why ATHAD stores both (README items 39, 42).

    max tau_layer = 0.254 at 3316 m       min 0.0099
    max tau_above = 4.40 at the surface   min 0.000 at the lid
    photosphere (tau_above = 1): 181 of 181 columns, median 3692 m, range 2637-4894 m

**No layer here carries `d(tau)` of order 1**, so the emission level is spread over ~4 layers
and the two-stream sweep — first order in `d(tau)` — resolves it. The same instrument in ATHAD
found **one cell carrying `d(tau)` = 55**, and later `tau_above` going 1626 -> 0.497 across a
single layer. That difference is a property of the atmospheres, not of the code: 250 bar of
water vapour puts the opaque-to-transparent transition where no reasonable grid resolves it.
Every column here has its photosphere inside the domain, against ATHAD's arms where **100 % of
columns radiate from the prescribed lid**.

### `div(u)` — and the projection does not achieve its target

**Read the call sites before the numbers.** In the TIME LOOP `PressureSolverAtm::run()` computes
`p_dyn` and **nothing else** — it never touches the velocity. The wind feels the pressure
through the `-dp/dr` term of the NEXT RK4 stage, and `pressure_stride = 4`
(`cAtmosphereModel.cpp:1020`) means even that happens on one iteration in four. The only place
an explicit correction `v <- v - grad(p)` is applied is `project_initial_velocity`, Step 3.

Measured there, in the model's own metric, term for term as the Poisson source builds it:

| `ATM_PROJ_SWEEPS` | div(u) rms | radial term | rms / radial |
|---|---|---|---|
| 1 (default) | 3.114e-03 | 3.633e-03 | **0.857** |
| 10 | 3.521e-03 | 4.090e-03 | **0.861** |
| 100 | 3.843e-03 | 4.393e-03 | **0.875** |

**The projection removes about 14 % of the divergence, and 100x the relaxation does not help.**
This is ATHAD item 72's finding — "converged to a fixed point that is not divergence-free" —
reproduced here with the direct instrument instead of inferred from `Psi`, and it explains why
`ATM_PROJ_SWEEPS` was inert on `Psi(ground)` above: the knob is not failing to reach a good
projection, the projection converges to something that is not divergence-free.

**Consequence for any vertical-velocity question.** The saved velocity field has been through a
full RK4 step plus the polar, orographic and radial filters since it was last projected, so its
divergence measures the filters rather than the solver, and **the vertical wind in the output is
not slaved to the horizontal flow**. A continuity test on saved output — the obvious way to ask
"is this vertical velocity right" — cannot answer the question here. That also accounts for a
3-5x mismatch found between the model's zonal-mean vertical wind and the wind its own
streamfunction implies, without needing the metric to explain it.

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
is a **null on every integrated quantity**: the residuum and the closure ratio both fall
slightly. Two trees, measured the same way, opposite outcomes — so *"the correct Jacobian makes
the projection worse"* is an ATHAD property and not a family law. This tree has the WORSE metric
(23.2x against 11.8x) and its INTEGRATED diagnostics do not notice the repair at all.

> ### CORRECTION, 2026-08-26, and it changes what the table above means
>
> The commit that landed this section (`5d46fb1`) read the -69 % radial-wind extremum as *"a
> spike being suppressed rather than a change in the circulation"*. **THAT IS WITHDRAWN. It was
> wrong, and the check that caught it was asked for by the reader of the write-up, not by its
> author.**
>
> The two arms were re-run with the full field dumped and the radial wind compared on two
> orthogonal slices (zonal, 41 levels x 181 latitudes; longitudinal, 41 x 361), ratios
> exact/legacy:
>
> | slice | RMS | p50 | p90 | p99 | max |
> |---|---|---|---|---|---|
> | zonal | **0.611** | **0.460** | 0.626 | 0.637 | 0.631 |
> | longitudinal | **0.537** | **0.597** | 0.506 | 0.493 | 0.491 |
>
> **The MEDIAN falls as far as the maximum, on both cuts, and every level between the
> boundaries falls by 0.50-0.69.** That is not a clipped extremum; the vertical-wind field is
> **40-50 % weaker everywhere**.
>
> So `ATM_METRIC_EXACT` is **not a null** — it is a change that the diagnostics used to judge it
> cannot see. `Psi` is built from the MERIDIONAL wind, KE is dominated by the horizontal
> components, and neither looks at the vertical wind, which is precisely the field a RADIAL
> metric governs. Same trap as ATHAD item 81's `Psi_max` false null, and the same lesson as its
> item 74: before believing a field did not move, check that the instrument reads it.
>
> **AND THE ATTRIBUTION QUESTION HAS NO ANSWER, BECAUSE IT PRESUMES A SPLIT THAT DOES NOT
> EXIST.** This block first named two candidates — (a) the exact Jacobian, `exp_rm` being 21.8x
> larger at the surface, and (b) the curvature term `-curv * df/dr` added with it — and called
> for one arm with the curvature forced off to separate them. `ATM_METRIC_NOCURV` was built and
> the arm was run. **Both the framing and the test design were wrong, and the run says so.**
>
> Difference RMS against the legacy arm, as a percentage of the legacy field, iteration 1:
>
> | | zonal | longitudinal |
> |---|---|---|
> | Jacobian only (`nocurv − legacy`) | 92.7 % | 83.1 % |
> | curvature only (`exact − nocurv`) | 72.8 % | 65.3 % |
> | **both together (`exact − legacy`)** | **60.8 %** | **51.8 %** |
>
> **The combined change is SMALLER than either piece alone, on both slices.** The two halves
> partially cancel, so there is no additive split to attribute — that is a statement about the
> arithmetic, not a hedge about the noise.
>
> **A second premise failed too, and it was mine.** The test was justified by "all three arms
> start from an identical state, so step 1 isolates the operator". They do not:
> `project_initial_velocity` uses `exp_rm` in the Poisson operator AND in the gradient
> correction, so **the arms already differ at iteration 0, before any time step** — by 95 % of
> the field, with RMS ratios 0.71 and 0.67. The metric acts through the INITIAL PROJECTION
> first, and the 40-iteration comparison above was never measuring divergence from a shared
> state either.
>
> **What the run does establish, and it is more useful than the attribution would have been:**
> the difference is **50–95 % of the field while the field's own RMS changes only 20–48 %**.
> Fields of comparable size, different SHAPE. `ATM_METRIC_EXACT` does not damp the vertical
> wind; it produces a **different vertical velocity field**, from iteration 0 onward, and the
> RMS halving is a downstream consequence rather than the effect. That is what the strange
> per-level profile was saying — 0.068 at one level and 1.20 at another is a RELOCATED
> structure, not a damped one.
>
> **This is therefore not decidable by more A/B runs.** It is a physics judgement — which field
> is right — and this tree has no instrument that can make it: `Psi` is built from the
> meridional wind, KE from the horizontal components, and the vertical wind is exactly what
> nothing integrates. Porting ATHAD's `div(rho u)/rho` print is the prerequisite for going
> further, not another arm.

**The default stays OFF, and the reason is no longer run length.** It is that the knob produces
a DIFFERENT vertical velocity field from iteration 0 onward, that no integrated diagnostic in
this tree registers it, and that the tree has no instrument able to say which field is right.
Forty iterations being a spin-up is a second reason and now the smaller one — ATHAD's item 54 watched a 4.4 % gap grow to 9.1 % between iterations
20 and 40.

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
| `ATM_CONV_ADJ` | `0` (off) | Dry convective adjustment, Manabe-Strickler, ported from ATHAD 2026-08-27. **This tree had none**, so nothing could remove a superadiabatic layer; see the section below. `ATM_CONV_ADJ_LAPSE` scales the critical lapse (`0` = isothermal criterion, `>1` stricter than dry) and `ATM_CONV_ADJ_PASSES` caps the sweeps per column (default 64) |
| `ATM_V_MASSBAL` | **`1` (ON since 2026-08-28)** | Removes the column-mean mass flux from the prescribed initial `v`, so `INT(rho*v*dz) = 0` per column. The prescribed cell is a linear ramp in height and closes in neither volume nor mass; see the streamfunction section. **`0` restores the old unbalanced profile exactly** — every measurement in this file dated before 2026-08-28 was made on that branch |
| `ATM_PROJECT_IN_LOOP` | `0` (off) | `<sweeps>`: a real velocity projection inside the time loop — seed `aux` from `u/v/w`, relax, apply `v <- v - grad(p)`, with `p_dyn` saved/restored. Distinct from the solver's own `run()`, which projects the momentum TENDENCY, not the velocity. Measured a null on `Psi(ground)` at 10 sweeps AND at 200 (-0.013 %) |
| `ATM_ANELASTIC` | `0` (off) | Solves `div(rho_bar u) = 0` instead of `div(u) = 0`. Null on a full run (-0.006 %) because the time loop never applies the pressure to the velocity; measured **2.0x the volume projection** across the INITIAL projection, where it is applied |
| `ATM_RAD_TOPO` | **`1` (ON since 2026-08-28)** | Puts the radiation column on `i_topography` instead of level 0. Over topography the sub-surface cells carry a real `p_stat`, so they enter `sum_dp` and dilute every air layer above them; and level 0, which `BC_Atm` stuffs with the mountain-top humidity while the rock stays dry, collects a large share of the water-vapour optical depth. Off-branch at 100 iterations: max `tau_layer` **7.885** over the Himalaya against 0.022 over ocean. `0` restores the sea-level column |
| `ATM_RHIE_CHOW` | `0` (off) | Fourth-difference pressure smoothing in the Poisson source, against the collocated-grid checkerboard. **Measured a null on `Psi(ground)` here (+0.005 %)** — it annihilates smooth fields by construction, so it cannot act on a domain-scale quantity |
| `ATM_CLOUD_TAU_MAX` | **`2.0` (ON since 2026-08-28)** | Per-layer ceiling on the CLOUD optical depth, scaling `LWP_i` and `IWP_i` together so the LW (`tau_cloud`) and SW (albedo bump) stay balanced. `cwp_cap_col` bounds the COLUMN condensate path and does not bound a LAYER: without this the shipped branch reaches layer emissivity **0.99962** aloft with 1773 cells above 0.9 per latitude slice — the near-blackbody pathology the de-saturation split exists to prevent. Clear-sky OLR is bit-identical across the flip (180.33882 W/m2), cloudy OLR +0.053. `0` disables |
| `ATM_PSI_PROJ_DUMP` | `0` (off) | Writes the streamfunction either side of `project_initial_velocity`, as iterations `-1` and `-2`. Print/CSV only |
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

### Dry convective adjustment (`ATM_CONV_ADJ`, 2026-08-27)

Until this port the tree had **no dry convective adjustment at all** — no such file, no call
site, only `MoistConvection.h`. Nothing could remove a superadiabatic layer, so one formed and
stayed.

Measured over 100 iterations at 24 threads. At iteration 0 the surface layer is *stable*
(-3.93 K/km, `brunt_N2 < 0` in 0.000 of columns). By iteration 20 the lowest 39 m runs at
**-18.05 K/km, twice the dry adiabat**, with `brunt_N2 < 0` in 59 % of columns at `i = 0` and
`i = 1` and 0 % at `i = 2` — one unstable layer, exactly the surface-to-first-level step. The
time loop makes it: over the run the surface warms +0.76 K while the air above warms +0.15 K,
and nothing mixes them.

With `ATM_CONV_ADJ=1`:

| quantity | off | on |
|---|---|---|
| surface lapse | -19.45 K/km | **-9.76** (dry adiabat -9.8) |
| `brunt_N2` median `i=0` | -2.69e-04 | **-1.0e-07** |
| `brunt_N2` median `i=1` | -4.75e-05 | **+6.09e-05** (frac<0 0.587 -> 0.006) |
| `i=3` and above | | unchanged to four figures |
| enthalpy drift | | **3.07e-16** |
| `Psi(ground)` / `max\|Psi\|` | | -0.02 % / -0.3 % |

**Measured again in isolation** (`ATM_CONV_ADJ=1` alone, `ATM_V_MASSBAL` off, against the
all-defaults reference), because the table above was measured on top of the mass balance and its
effect on `Psi` was therefore confounded with a field that repair had already changed by 71 %:

| quantity | reference | conv-adj alone |
|---|---|---|
| `brunt_N2` median `i=0` | -2.6811e-04 | **-1.0e-07** |
| `brunt_N2` median `i=1` | -4.706e-05 | **+6.091e-05** |
| `brunt_N2` frac<0 `i=1` | 0.58726 | **0.00554** |
| lapse 0-1 | -19.362 K/km | **-9.758** |
| `Psi(ground)` rms | 1.6657e+11 | 1.6634e+11 (**-0.139 %**) |
| `max\|Psi\|` | 4.1568e+11 | 4.1532e+11 (-0.087 %) |

**The thermal result is identical to the on-top-of-massbal arm** — `brunt_N2` at `i=0` lands on
-1e-07 in both, the lapse on -9.758 against -9.76, `frac<0` at `i=1` on 0.0055 in both — and the
scheme statistics differ by two columns and seven layers out of ninety-six thousand. So the
adjustment's effect does not depend on the initial mass balance.

**The two repairs are orthogonal, and the obvious worry does not materialise.** On the
uncorrected profile the bottom layers start ~10 K/km more unstable, so mixing them changes
density more and might have been expected to move `Psi`; it moves it **-0.139 %**, seven times
the -0.02 % measured on the corrected field and still negligible. `Psi_max` stays at 15N, z = 0 m
at all five printouts either way — the convective adjustment does nothing about the ground
non-closure, which is `ATM_V_MASSBAL`'s job.

67.3 % of columns adjust, worst column 2 sweeps of an allowed 64. It runs **before**
`densities()`, which is where `brunt_N2` is formed — otherwise the instrument reports the profile
the adjustment was about to remove.

**Two numbers here must be read by magnitude, not by count.** `frac<0` at `i = 0` stays 0.590
while the value falls to -1e-07: that is the scheme's tolerance on a layer now neutral by
construction, so the fraction counts round-off there. And max/min temperature are *bit-identical*
between the arms while the median surface profile moves 0.3 K — the extremum is a false null.

The port kept ATHAD's per-layer critical drop from `get_layer_height()` and its cumulative-sum
segment adiabat, because the family's shared `planet/ConvectiveAdjustment.h` computes
`dz_m = L_atm*1e3/(im-1)`, a layer thickness only on a **uniform** grid. Two things were adapted:
the column starts at `i_topography` rather than level 0 (ATHAD has no terrain), and `cp` is the
constant `cp_l` rather than a mixture value (this tree has no `MixtureAtm`).

### Meridional streamfunction as a field

`Psi` — the meridional mass streamfunction [kg/s] — is published as a 3D field as well as a
printed extremum. `write_meridional_streamfunction()` already formed `psi[i][j]` and threw it
away except for a CSV and the `Psi_max` line; it now fills `Psi`, which appears in
`print_min_max_atm` and in the radial/zonal/longal ParaView dumps as `PsiMerid`. It is a zonal
mean, so it is replicated across `k`.

**`Psi` does not close at the ground, and most of it is the prescribed initial profile.**
`Psi(lid)` is `0.0000e+00` exactly at every iteration while `Psi(ground)` rms is 1.55e+11 kg/s —
about 40 % of the Hadley cell's own strength — and `u` is pinned to zero at BOTH walls
(`BC_Atm.h:679`, `:697`, verified in the field, `max|u| = 0.000000e+00`). With no flux through
either wall, `div(rho v) = 0` would force `Psi(ground) = 0`. It is not a boundary cell: rms`|Psi|`
decays smoothly through the whole depth. 93.8 % of it is present at iteration 1, and
`project_initial_velocity` changes it by **-0.014 %**, so the projection does not remove it.

The cause is `VelocityInitializer::init_v_or_w`, which builds `v` as a **linear ramp in height**
between two hand-set endpoints. Nothing constrains `INT(rho*v*dz) = 0`, or even `INT(v*dz) = 0`:
at 15N the column runs +3.67 m/s at the ground against -0.54 m/s at 7.4 km, and `INT(v*dz)` is
+9151 m^2/s, so it fails in volume as well as mass. `ATM_V_MASSBAL=1` subtracts the
density-weighted column mean, which zeroes the column mass flux while shifting the profile by a
single constant, leaving the shear — and hence the overturning — untouched. **-94.8 % at
initialisation, -71.2 % at iteration 100**, with `max|Psi|` 4.16e+11 -> 1.61e+11. It is **not a
cure**: the residual rises monotonically and no limit is claimed.

The reason for publishing it is ATHAD's README item 68: there `Psi_max` was reporting the
**spurious surface mass flux** rather than the circulation, which a scalar cannot show and a
field can. **The fill now runs BEFORE `print_min_max_atm`** — the other order left the reported
`Psi` min/max one checkpoint stale and zero on the first, the same defect ATHAD's item 42
records. The CSV and the vtk were always correct; only the printed extrema were behind.

### Pressure-field diagnostics

`reportDivergence()` prints, at both the initial projection and each pressure solve: the
divergence the projection is meant to remove, a **checkerboard index** on `p_dyn`
(rms(`p` - mean of the 6 neighbours)/rms(`p`); 0 smooth, 2 pure Nyquist), the **Nyquist share of
the anomaly along each axis** with absolutes beside every ratio, and the **Poisson operator
weights** `c_r`/`c_the`/`c_phi`.

Read the per-axis absolutes, not the global index. The global index divides by rms(`p_dyn`) and
so cannot be compared between trees, between arms, **or between iterations of the same run** —
in ATHAD the same field read 0.0003 globally and 0.961 on the zonal axis, and here the index
falls 0.63 -> 0.0384 over 100 iterations purely because its denominator grows 24x, while the
absolute grid-scale amplitude *rises* 2.99e-06 -> 4.45e-06. A falling index is not a decaying
mode.

Neither instrument isolates the Nyquist mode from smooth curvature — both are second differences
— so a growing smooth field contaminates both. Of the two, the per-axis *share* is the most
trustworthy.

### ParaView field dumps

The atmosphere's ParaView writers (`atmosphere/Paraview_Atm.cpp`) expose many optional fields
as commented-out `dump_*` lines — force components (`BuoyancyForce`, `CoriolisForce`,
`CentrifugalForce`, `PresGradForce`), `CO2-Concentration`, latent/sensible/radiative heat, and
more. Uncomment the ones you need for a given investigation and rebuild. The pole-singular
`paraview_sphere_vts` output can be skipped in favour of the `panorama` `.vts` (whose
derivatives are well-behaved away from the poles).

## 2026-08-28: the projection arms closed out, and the surface boundary condition

### `ATM_PROJECT_IN_LOOP` at 200 sweeps — still a null, and now the null means something

`f88d319` recorded the knob as implemented-but-untested and discounted its 10-sweep null on the
grounds that the projection was under-resourced: `project_initial_velocity` runs 200 relaxation
passes and the in-loop call had been given 10, with red-black Jacobi converging the largest
scales slowest — error decaying like `(1 - c/N^2)` per sweep, and `Psi(ground)` exactly such a
domain-scale mode. It has now had 200. Iteration 100, 24 threads, `nm = 100`:

| sweeps | `Psi(ground)` rms | vs off | `div(u)` at solve 25 | wall |
|---|---|---|---|---|
| off | 1.66571090e+11 | — | 1.766e-03 | 527 s |
| 10 | 1.66569245e+11 | -0.0011 % | 1.764e-03 | 556 s |
| 10 + anelastic | 1.66561990e+11 | -0.0055 % | — | 552 s |
| **200** | **1.66549029e+11** | **-0.0132 %** | **1.761e-03** | **1017 s** |

**Twenty times the sweeps bought 2.5x the `div(u)` reduction** (-0.28 % against -0.11 %) and left
`Psi(ground)` where it was. The cost estimate in the first write-up was wrong too: "~20x the
solver cost per iteration" is true per solve, but the solve is about a fifth of the runtime —
**1.93x wall clock, 17 minutes**. The arm was deferred as expensive and was not.

**What it does not settle, and the gap is one print.** `div(u)` is only ever reported AT THE
PRESSURE SOLVE, never immediately after `project_velocity_in_loop`. "The projection zeroes
`div(u)` and one time step puts all of it back" and "the projection never reduces `div(u)`" have
the identical signature in that table, and they are different defects.

### `Psi(lid) = 0` is the integration constant, and `im` cannot close `Psi(ground)`

`MinMax_Atm.cpp:232` integrates DOWNWARD from the lid with `psi[im-1] = 0` by construction, so
the `Psi(lid)` column of the streamfunction table above is arithmetic and is identical in every
arm ever run. The earlier framing, "the lid closes exactly and the ground does not", gave it
weight it never had. **All the content is in `Psi(ground)`** — which is therefore identically the
whole-column integral `2*pi*a*cos(phi)*INT(rho*vbar dz)`, not a value AT `i = 0` but a TOTAL
reported at `i = 0`.

That reading rules out two tempting non-cures:

1. **Reshaping the initial `v(z)` to push the offset up the column cannot work.** A definite
   integral is invariant to any vertical redistribution preserving `INT(rho*v dz)`: the profile
   can be reshaped until `Psi` is near-zero through most of the depth and the ground value will
   not move by one digit. Re-integrating upward from the ground is the cosmetic version — it
   relabels which end carries the constant and destroys the instrument.
2. **More vertical levels cannot close it, because the defect is in the INTEGRAND.** Measured by
   re-integrating the written field with a cubic spline on the SAME 41 levels (nodal `rvbar`
   recovered from the `psi` differences; that back-out is marginally stable, so the lid anchor
   was swept over its plausible range and the three anchors agree to three decimals):

   | field | trapezoid (model) | cubic spline | change |
   |---|---|---|---|
   | shipped, iter 20 | 1.622498e+11 | 1.619763e+11 | **-0.169 %** |
   | shipped, iter 100 | 1.665711e+11 | 1.663176e+11 | **-0.152 %** |
   | `ATM_PROJECT_IN_LOOP=200`, iter 100 | 1.665490e+11 | 1.662955e+11 | -0.152 % |
   | `ATM_V_MASSBAL=1`, iter 100 | 4.793468e+10 | 4.777011e+10 | -0.343 % |

   Quadrature is ~0.15 % of a 40 % non-closure, and in absolute terms ~1.6e8-2.5e8 kg/s in every
   arm — so even the residual `ATM_V_MASSBAL` leaves is **200x larger than the discretisation
   error**. The source says the same without running anything: `VelocityInitializer.h:55` sets
   the 15N Hadley column with `init_v_or_w(m.v, 75, -3.0, 4.0)`, tropopause **-3.0** and surface
   **+4.0**, and for the linear ramp the continuous integral is `H*(v_s+v_t)/2 ∝ +1.0`, nonzero
   before any grid exists.

   *Caveat kept rather than argued away*: this tree's stretch is the family's worst (`zeta`
   3.715, spread 23.21x, 38.9 m bottom layer), so a real trapezoid error could have hidden near
   the surface. The spline measurement is what bounds it at 0.15 %; it is a bound, not an
   argument.

### `ATM_V_MASSBAL` is the default since 2026-08-28

No other lever in this tree touches `Psi(ground)`: `ATM_PROJ_SWEEPS` is inert here, `ATM_RHIE_CHOW`
is structurally incapable of it, `ATM_ANELASTIC` cannot reach the velocity in the loop, and
`ATM_PROJECT_IN_LOOP` is a null at 10 and at 200 sweeps. The column mass balance removes 94.8 % at
initialisation and 71.2 % at iteration 100. Verified both directions at iteration 100, `nm = 100`,
24 threads:

| branch | `Psi(ground)` | `max\|Psi\|` | matches |
|---|---|---|---|
| `ATM_V_MASSBAL=0` | 1.6657e+11 | 4.1568e+11 | the recorded old default |
| default (unset) | 4.7935e+10 | 1.6078e+11 | the recorded `=1` arm |

**And the global maximum of `Psi` moves off the ground**, which is ATHAD item 68's signature:

| | max signed `Psi` | min signed `Psi` | global `max\|Psi\|` |
|---|---|---|---|
| massbal OFF | +4.1568e+11 @ 15N, **0 m** | -3.5102e+11 @ 45S, 1729 m | 4.1568e+11 @ 15N, **0 m** |
| massbal ON | +1.4689e+11 @ 15N, 0 m | -1.6078e+11 @ 45S, 3316 m | **1.6078e+11 @ 45S, 3316 m** |

On the old default the largest value in the streamfunction WAS the ground defect; on the new one
the real circulation exceeds it. **Every other knob measurement in this file dated before
2026-08-28 was made against the unbalanced initial state**, so those comparisons now differ in
their initial state as well as in the quantity under test — the trap this project already records
ATHAD losing an attribution to. `ATM_V_MASSBAL=0` restores the old branch exactly.

### `brunt_N2 < 0` at `i = 0` is the OCEAN MASK, and `t` has no radial boundary condition there

Measured on the level-0 radial slice at iteration 100:

| | cells | `P(N2 < 0)` | median `N2` |
|---|---|---|---|
| ocean | 43 602 (66.73 %) | **0.9998** | **-3.70e-04** |
| land | 21 739 (33.27 %) | **0.0049** | **+2.11e-04** |

43 593 of the 43 700 unstable cells are ocean and 107 are land. **The "59 % of columns" recorded
for `ATM_CONV_ADJ` was always the ocean fraction.**

Four facts, and the last is the defect:

1. **Level 0 is not prognostic** — `RungeKutta_Atm_Turb.cpp:111` is `for(int i = 1; i < im-1;)`.
2. **`t` has NO radial BC at `i = 0`.** `bcRadius` carries three pattern lists — `both_cubic`,
   `vn_bot_cubic_top`, `cubic_bot_vn_top` — and `t` is in **none** of them, while `p_stat`,
   `r_humid`, `cloud`, `ice` and `tke` all get one. `t` has explicit special handling at the LID
   (`pin_t_top`) and nothing at the bottom.
3. **The mechanism meant to cover that is a self-assignment over ocean.** `BC_Atm.h:433`,
   documented as "copy surface-level (`i_mount`) values down to the i=0 reference layer … so
   surface fluxes see the correct surface state", is `t.x[0] = t.x[i_topography]`. Over land
   `i_topography > 0` and it does real work — which is exactly why land is stable. Over ocean
   `i_topography == 0` (`FileIO_Atm.cpp:488`) and it reads `t.x[0] = t.x[0]`.
4. **So level 0's only coupling to the air above is a Shapiro smoother**, `damp_wiggles(t,…)` at
   `cAtmosphereModel.cpp:1113`. Level 0 drifts **+0.786 K** over 100 iterations against +0.15 K
   above it.

**The stencil explains why it is ONE layer.** `N2` at `i = 0` is one-sided on `theta[0],
theta[1]`; at `i = 1` centred on `theta[0], theta[2]`; at `i = 2` centred on `theta[1],
theta[3]` — the first stencil that does not see level 0. The 59/59/0 pattern is one bad value
read by two stencils, not two unstable layers.

**The model has not decided whether ocean level 0 is the ocean or the air.** If it is a
prescribed SST skin then `N2 < 0` there is physically correct and the repair is to credit the
bulk flux to level 1; if it is the lowest air level it needs RK4 and a flux BC. It is currently
neither. **This reframes `ATM_CONV_ADJ`**: its -19.45 -> -9.76 K/km is a palliative mixing away a
boundary-condition gap, not a missing convection scheme.

**And `Q_Sensible` is NOT the surface flux** — an earlier note ran two quantities together. The
array at `RHS_Atm_Turb.cpp:513` is `coeff_S*lap(T)`, the CONDUCTIVE flux divergence at every
interior cell; crediting it anywhere would double-count the thermal diffusion `diffusion_t_re`
already applies. The surface bulk flux is the separate `c_H*(T_s - T_air1)` with
`c_H = 15 W/m2/K` at `MultiLayerRadiation.h:425`. Both are real defects; only the second bears on
the instability, and MLR's copy is refreshed only every 20 iterations. **An earlier version of
this section said MLR never runs inside the loop, citing "calls at `cAtmosphereModel.cpp:713-889`,
loop at 985". That was wrong**: 711-892 are LAMBDA DEFINITIONS and the call sites are at 1463-1471,
inside the loop, with the active `radiation_mode` 5 invoking `cloud_radiation_diag()` — and hence
MLR — every `teq_refresh_stride` = 20 iterations.

`ATM_SFC_FLUX=<c_H>` (default `0` = off, bit-identical; `15` matches MLR) forms that flux in the
loop from the live `t`, as a rate `k_S = c_H/(rho*cp*dz)` non-dimensionalised in advective time
`k_S*L_atm/u_0` exactly like the Held-Suarez relaxation, applied at the first air level
`i_topography+1`. The surface is treated as a fixed reservoir — the prescribed-SST convention,
stated rather than hidden; crediting the surface back requires level 0 to be prognostic, which is
the larger question this defect raises.

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
