# OxySim-129

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23240592.svg)](https://doi.org/10.5281/zenodo.23240592)

**A High-Fidelity CFD Framework for Carbonaceous Solid Fuel Combustion**
*From Single-Particle Studies to MW-scale Reactor Simulations*

OxySim-129 is an OpenFOAM-based C++ framework for Lagrangian particle simulations in reactive flows with flamelet-based tabulated chemistry. It provides solvers and libraries for coal and biomass combustion using the carbonaceous particle cloud coupled to flamelet manifolds, covering a range of scales from small particle groups to MW-scale swirl-stabilized combustors. Thermochemistry is looked up from tabulated flamelet manifolds via [`flut-reader`](https://github.com/STFS-TUDa/flut-reader), STFS's companion FLUT library (included as a submodule, see [Installation](#installation)).

Developed at the [Institute for Simulation of reactive Thermo-Fluid Systems (STFS)](https://www.stfs.tu-darmstadt.de), TU Darmstadt, within the [DFG CRC/Transregio 129 "Oxyflame"](https://www.oxyflame.de/) (project number 215035359).

## Contents

- **[For users](#for-users)**
- **[Installation](#installation)**
- **[Running](#running)**
- **[Solvers](#solvers)**
- **[Code structure](#code-structure)**
- **[Testing](#testing)**
- **[References](#references)**

### Component documentation

| Component | Description |
|---|---|
| [flameletCloudFoam](applications/solvers/flameletCloudFoam/README.md) | Flamelet + Lagrangian particle solver |
| [flameletEquations](src/flamelet/flameletEquations/README.md) | Scalar transport equations + EquationHandler |
| [Thermo](src/thermo/basicThermoSTFS/README.md) | Base thermo class, fans out to cantera and flamelet thermo |
| [flut-reader](https://github.com/STFS-TUDa/flut-reader) | Flamelet look-up table (FLUT) library (submodule) |
| [CloudHandler](src/CloudHandler/README.md) | Runtime cloud manager |
| [CarbonaceousParcel](src/lagrangian/intermediate/parcels/Templates/CarbonaceousParcel/README.md) | Carbonaceous parcel template |
| [Submodels](src/lagrangian/intermediate/submodels/README.md) | Devolatilization, heat capacity, char oxidation, injection, drag, cloud radiation |
| [Radiation](src/radiation/README.md) | WSGG / fvDOM radiation model |
| [test\_fp](tests/integration/test_fp/README.md) | Freely-propagating flame integration test |
| [test\_C2SM](tests/integration/test_C2SM/README.md) | C2SM devolatilization integration test |
| [test\_flamelet\_particle\_heating](tests/integration/test_flamelet_particle_heating/README.md) | Particle heat transfer integration test |
| [test\_flamelet\_particle\_gas\_phase\_coupling](tests/integration/test_flamelet_particle_gas_phase_coupling/README.md) | Flamelet–particle coupling integration test |
| [test\_multiphaseRadiation\_dorigon](tests/integration/test_multiphaseRadiation_dorigon/README.md) | Gas+particle radiation integration test (air atmosphere) |
| [test\_multiphaseRadiation\_bordbar](tests/integration/test_multiphaseRadiation_bordbar/README.md) | Gas+particle radiation integration test (oxyfuel atmosphere) |

---

## For users

OxySim-129 is maintained by a subset of our group alongside our regular research work. Like OpenFOAM itself it comes with no warranty of any kind. If something looks like a bug, we'd genuinely like to hear about it: please open a ticket on the ["Issues"](https://github.com/STFS-TUDa/oxysim-129/issues) page.

We try to answer, but given our team size we can't provide dedicated support for every setup question.

For anything you'd rather discuss directly with the authors, reach us at steffens@stfs.tu-darmstadt.de or vahl@stfs.tu-darmstadt.de.

---

## Installation

### Dependencies

| Dependency | Version |
|---|---|
| [OpenFOAM](https://www.openfoam.com/) | [v2512+](https://gitlab.com/openfoam/core/openfoam/-/tree/OpenFOAM-v2512?ref_type=tags) |
| Cantera | 3.1+ |
| HDF5 | C library + headers, required by [flut-reader](https://github.com/STFS-TUDa/flut-reader) (see its README) — commonly already provided alongside OpenFOAM, since OpenFOAM itself links HDF5 |
| CMake | 3.28+ |
| C++ compiler | GCC 13 or equivalent (C++20) |

### Clone

```bash
git clone --recurse-submodules https://github.com/STFS-TUDa/oxysim-129.git
cd oxysim-129/
```

The `--recurse-submodules` flag is required — it fetches the flamelet look-up table library (`extern/flut-reader`) that the solvers depend on.

Git LFS (https://git-lfs.com/) is also required (git lfs install, once per machine) — reference/test data (*.h5, *.png, *.parq, mesh files, etc.) is stored via LFS, and without it you won't get the real files.

If you already cloned without it, initialize the submodule afterwards instead of re-cloning:

```bash
git submodule update --init --recursive
```

### Set up the Python environment

The Python environment is required for running the integration tests. A virtual environment keeps its packages (`pytest`, `numpy`, `pandas`, `scipy`, `matplotlib`, `pyvista`) out of your system/global Python install:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r tests/python_utils/requirements.txt
```

Any other conda/venv/module environment providing those same packages works just as well — there is no dedicated install script beyond `requirements.txt`.

### Build

**1. Load the OpenFOAM environment**

OpenFOAM's environment must be active before running CMake, along with Cantera and CMake at the versions listed above.

- If your system provides these via a module system (e.g. Environment Modules, Lmod — common on HPC clusters), load them the way your site names them, e.g. `module load <your-openfoam-module> <your-cantera-module> <your-cmake-module>`. Module names are site-specific, so check what's available with `module avail`.
- Otherwise, source OpenFOAM's environment directly and point CMake at your own Cantera installation (see step 2):
  ```bash
  source /path/to/OpenFOAM-v2512/etc/bashrc
  ```

For convenience, [`scripts/source_template.sh`](scripts/source_template.sh) bundles this together with the `CANTERA_ROOT` and [Running](#running) environment setup below into a single script — copy it to `scripts/source_local.sh` (gitignored), fill in your paths, and `source` it instead of repeating these steps by hand.

**2. Configure and build**

```bash
mkdir -p build && cd build
cmake .. -DCANTERA_ROOT=/path/to/cantera
make -j$(nproc)
```

If a module system set `CANTERA_ROOT` for you in step 1, the `-D` flag above can be omitted. Otherwise, replace the path with Cantera's installation prefix — the root directory that contains `include/` and `lib/` subdirectories.

**Finding the path on your system:**
- If installed via a module system: loading the module typically sets `CANTERA_ROOT` automatically, and no `-D` flag is needed.
- If installed via a package manager (apt, brew, Spack, conda): run `find /usr /opt $HOME -name "cantera.pc" 2>/dev/null` to locate the Cantera prefix, then strip `/lib/pkgconfig/cantera.pc` from the result.
- If built from source: use the directory passed as `--prefix` or `CMAKE_INSTALL_PREFIX` during that build.

The build output goes to `build/plattforms-stfs-foam/<platform>/default/`:
- Solvers → `bin/` (`flameletFoam`, `flameletCloudFoam`)
- Libraries → `lib/`

#### Build flags

| CMake variable | Effect |
|---|---|
| `COMPILE_LAGRANGIAN_WITH_TESTS=TRUE` | Compile Catch2 unit tests |
| `OXYSIM_INSTALL_DIR=<name>` | Change the output sub-directory (default: `default`) |

---

## Running

Building only produces the binaries — running a case additionally requires `PATH` and `LD_LIBRARY_PATH` to point at the build output, so the solvers are found and can dynamically link against their libraries.

With the OpenFOAM environment still loaded (see [Build](#build)):

```bash
export FOAM_USER_APPBIN="$(pwd)/build/plattforms-stfs-foam/${WM_OPTIONS}/default/bin"
export FOAM_USER_LIBBIN="$(pwd)/build/plattforms-stfs-foam/${WM_OPTIONS}/default/lib"
export PATH="${FOAM_USER_APPBIN}:${PATH}"
export LD_LIBRARY_PATH="${FOAM_USER_LIBBIN}:${LD_LIBRARY_PATH}"
```

Replace `default` with the value passed to `OXYSIM_INSTALL_DIR` at build time, if any.

If Cantera was loaded via a module, its libraries are already on `LD_LIBRARY_PATH`. If it was installed manually, add it too:

```bash
export LD_LIBRARY_PATH="${CANTERA_ROOT}/lib:${LD_LIBRARY_PATH}"
```

Verify the solver resolves to your build:

```bash
which flameletCloudFoam
```

Then run a case as usual, e.g. `flameletCloudFoam` from the case directory.

### Demo cases

[`oxysim-demo-cases`](https://github.com/STFS-TUDa/oxysim-demo-cases) has two runnable `flameletCloudFoam` configurations to try once the build works — real cases rather than synthetic tutorials: a laminar coal case with volatile release only, based on Nicolai et al. (2021) [[11]](#references), and a turbulent (LES) oxyfuel biomass case with volatile release, char oxidation, and methane-assistance, based on Vahl et al. (2026) [[12]](#references). See that repository's README for prerequisites (Git LFS) and per-case run instructions.

---

## Solvers

### `flameletFoam`

Pressure-based solver for reactive flows using flamelet-tabulated chemistry. Solves transport equations for the flamelet input variables (mixture fraction Z, progress variable Y_c, enthalpy h_a, mixture fraction variance Z'^2) and looks up thermochemical quantities from a FLUT file.

**Required case files:**
- `constant/flameletProperties` — scalar transport equations, table lookup settings
- A FLUT file (flamelet look-up table in HDF5 format)

### `flameletCloudFoam`

Extends `flameletFoam` with Lagrangian carbonaceous particle tracking. Particle-to-gas source terms (mass, momentum, energy, species) are added to all governing equations. Uses `CanteraThermo` to interface gas-phase species with the particle models.

**Required case files (in addition to `flameletFoam`):**
- `constant/<cloudName>Properties` — cloud type, injection, submodel choices
- `system/controlDict` → `clouds` dict:
  ```
  clouds
  {
      coalCloud    carbonaceousCloud;
      biomassCloud carbonaceousCloud;
  }
  ```
- A mechanism file listing relevant species (O2, volatile species, Z, etc.)
- Optional: `constant/deleteParticlesDict` for particle removal criteria

---

## Code structure

```
src/
├── flamelet/
│   ├── flameletThermo/           Flamelet thermophysical model (table lookup, field updates)
│   ├── flameletBCs/              Boundary condition for flamelet scalar fields
│   └── flameletEquations/        Scalar transport equations + EquationHandler (incl. particle coupling)
├── thermo/
│   ├── basicThermoSTFS/          Base thermo class (STFS extensions)
│   └── canteraThermo/            Cantera-backed species mixture for particle–gas coupling
├── benchmarkUtility/             Timing utility for solver performance profiling
├── CloudHandler/                 Runtime-selectable multi-cloud manager
├── radiation/                    WSGG / fvDOM radiation model (see below)
└── lagrangian/
    └── intermediate/
        ├── clouds/               CarbonaceousCloud template + base class + derived
        ├── parcels/              CarbonaceousParcel template + derived + submodel macros
        └── submodels/
            ├── Carbonaceous/     Devolatilization (C2SM), heat capacity, char oxidation
            └── CloudRadiation/   Particle-cloud WSGG radiation (wsggGreyCloudAbsorptionEmission)

applications/
└── solvers/
    ├── flameletFoam/
    └── flameletCloudFoam/

extern/
└── flut-reader/                  Flamelet Look-Up Table library (submodule)

tests/
├── integration/
│   ├── test_fp/                              Freely propagating flame (flameletFoam)
│   ├── test_standard_BC_lookup/              Flamelet boundary condition lookup
│   ├── test_C2SM/                            C2SM devolatilization validation
│   ├── test_flamelet_particle_heating/       Particle heat transfer
│   └── test_flamelet_particle_gas_phase_coupling/  Full flamelet–particle coupling
├── conftest.py                   Shared pytest fixtures (testdir, write_report)
└── python_utils/
    ├── foam_utils.py              OpenFOAM-to-pandas helper (foam_to_df) used by test_fp
    └── requirements.txt           Python environment for running the integration tests
```

### Carbonaceous submodels

**Devolatilization** (`submodels/Carbonaceous/DevolatilizationModelSTFS/`):
- `C2SM` — Kobayashi competing two-step model [[1]](#references). Kinetic parameters from pkp-lite [[4]](#references).
- `NoDevolatilisationSTFS` — no-op

**Heat capacity** (`submodels/Carbonaceous/HeatCapacityModelSTFS/`):
- `ConstantHeatCapacity` — fixed Cp
- `HeatCapacityFit` — polynomial: cp(T) = a + bT + cT² + dT³ + eT⁴

**Surface reactions** (`submodels/Carbonaceous/SurfaceReactionModelSTFS/`):
- `COxidationKineticDiffusionLimitedRateSTFS` — Baum–Street char oxidation [[2]](#references), coefficients from [[3]](#references)

### CloudHandler

Manages one or more `carbonaceousCloud` instances at runtime. Reads the `clouds` dict from `system/controlDict`, evolves all clouds each time step, and accumulates particle source terms for the solver equations.

Supports conditional particle deletion via `constant/deleteParticlesDict` (criteria: position x/y/z, radial distance rz, velocity Umag/Uz, temperature T).

### Radiation

`fvDOMwsgg` solves the radiative transfer equation with a non-grey weighted-sum-of-gray-gases (WSGG) model, selectable for oxyfuel (`wsggBordbarAbsorptionEmission`) or air (`wsggDorigonAbsorptionEmission`) atmospheres. `flameletCloudFoam` calls it once per PIMPLE iteration and couples the resulting source term into the gas-enthalpy equation; combined with `wsggGreyCloudAbsorptionEmission`, it also accounts for emission/absorption by Lagrangian carbonaceous particle clouds. See [Radiation](src/radiation/README.md) for configuration details.

---

## Testing

### Integration Tests

Run the integration tests after loading the environment (see [Build](#build)) and activating the Python environment (see [Set up the Python environment](#set-up-the-python-environment)):

```bash
# Flamelet solver tests
python -m pytest tests/integration/test_fp/ -v
python -m pytest tests/integration/test_standard_BC_lookup/ -v

# Flamelet–particle coupling tests
python -m pytest tests/integration/test_C2SM/ -v
python -m pytest tests/integration/test_flamelet_particle_heating/ -v
python -m pytest tests/integration/test_flamelet_particle_gas_phase_coupling/ -v

# Radiation tests
python -m pytest tests/integration/test_multiphaseRadiation_dorigon/ -v
python -m pytest tests/integration/test_multiphaseRadiation_bordbar/ -v
```

Each test copies `foam_template/` into a temporary directory, runs the OpenFOAM case, and compares results against stored reference data.

### Unit tests

The `flameletThermo` library has a C++ unit test suite (Catch2, fetched automatically by CMake), built separately from the solvers and run via its own `Allrun` script:

```bash
COMPILE_LAGRANGIAN_WITH_TESTS=TRUE cmake .. && make -j10 testFlameletThermo
export OXYSIM_REPO_ROOT=/path/to/oxysim-129
bash tests/unit/testFlameletThermo/case/Allrun
```

---

## References

1. Kobayashi et al. (1977) doi:[10.1016/S0082-0784(77)80341-X](https://doi.org/10.1016/S0082-0784(77)80341-X) — C2SM devolatilization model
2. Baum & Street (1971) doi:[10.1080/00102207108952290](https://doi.org/10.1080/00102207108952290) — char oxidation model
3. Gövert et al. (2017) doi:[10.1016/j.fuel.2017.03.009](https://doi.org/10.1016/j.fuel.2017.03.009) — char oxidation coefficients
4. Vascellari et al. (2013) doi:[10.1016/j.fuel.2013.06.014](https://doi.org/10.1016/j.fuel.2013.06.014) — pkp-lite kinetics approach
5. Steffens et al. (2024) doi:[10.1016/j.fuel.2024.132098](https://doi.org/10.1016/j.fuel.2024.132098) — coal/biomass flamelet model
6. Nicolai et al. (2021) doi:[10.1016/j.combustflame.2021.111722](https://doi.org/10.1016/j.combustflame.2021.111722) — laminar flow reactor validation
7. Berkel et al. (2025) — optically accessible solid fuel combustor validation
8. Trivic doi:[10.1016/j.ijheatmasstransfer.2003.09.027](https://doi.org/10.1016/j.ijheatmasstransfer.2003.09.027) — non-grey WSGG–discrete-ordinates radiation formulation
9. Bordbar et al. (2014) doi:[10.1016/j.combustflame.2014.03.013](https://doi.org/10.1016/j.combustflame.2014.03.013) — WSGG coefficients for oxyfuel atmospheres
10. Dorigon et al. (2013) doi:[10.1016/j.ijheatmasstransfer.2013.05.010](https://doi.org/10.1016/j.ijheatmasstransfer.2013.05.010) — WSGG coefficients for air atmospheres
11. Nicolai et al. (2021) doi:[10.1016/j.proci.2020.06.081](https://doi.org/10.1016/j.proci.2020.06.081) — pulverized coal particle group combustion; basis for [`oxysim-demo-cases`](https://github.com/STFS-TUDa/oxysim-demo-cases)' `laminar-coal` case
12. Vahl et al. (2026) doi:[10.1016/j.fuel.2026.138602](https://doi.org/10.1016/j.fuel.2026.138602) — LES of swirl-stabilised gas-assisted oxy-fuel biomass flames; basis for [`oxysim-demo-cases`](https://github.com/STFS-TUDa/oxysim-demo-cases)' `turbulent-oxyfuel-biomass` case

---

## Acknowledgements

This work has been funded by the German Research Foundation (DFG) – project number 215035359 – within the framework of the CRC/Transregio 129 "Oxyflame".

## How to cite

If you use this software, please cite it using the metadata in [`CITATION.cff`](CITATION.cff), or click "Cite this repository" in the GitHub sidebar. A DOI for citing a specific version is also available: [10.5281/zenodo.23240592](https://doi.org/10.5281/zenodo.23240592).
