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
