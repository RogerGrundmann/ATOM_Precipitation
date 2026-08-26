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
