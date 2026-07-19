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
| `--nm` | `400` | Iterations per atmosphere/hydrosphere run |
| `--checkpoint` | `100` | VTK/panorama printout stride, and the convergence-monitor sampling cadence |
| `--tol` | `0.05` | Stop when `max|ΔSST|` between rounds drops below this (K) |

Each round runs from scratch in its own sub-directory; the report-only convergence
monitor writes a `convergence.csv` (ocean-mean-T and kinetic-energy drift) per round.

## Output

Each run writes output files to `output/`:

- **VTK files** — radial, zonal, and longitudinal slices readable by [ParaView](https://www.paraview.org/)
- **XYZ grid files** — tab-separated gridded data for post-processing

## Python interface

Install the Python bindings:

```bash
pip install -e python/
```

Then use the `Atmosphere` and `Hydrosphere` classes:

```python
from atom import Atmosphere

atm = Atmosphere("cli/config_atm.xml")
atm.run()
```

See `benchmark/run.py` for a fuller example and `benchmark/Demo.ipynb` for visualisation.

## Input data

The `data/` directory ships with:

- Paleotopography/bathymetry grids 0–140 Ma (Smith et al. 1994; Golonka et al. 1997)
- Present-day surface temperature (NASA)
- Present-day precipitation (NASA)
- Present-day salinity (NASA)

## Authors

Roger Grundmann, Michael Chin
