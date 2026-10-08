# Basic Thermo STFS
↑ [OxySim-129](./../../../README.md)

[basicThermoSTFS.H](./basicThermoSTFS.H)
[basicThermoSTFS.C](./basicThermoSTFS.C)

## Overview

`basicThermoSTFS` is a thin base on top of OpenFOAM's `fluidThermo`, shared by [`flameletThermo`](./../../flamelet/flameletThermo) and [`canteraThermo`](./../canteraThermo). It gives every `fluidThermo`/`basicThermo` virtual a default `NotImplemented` body, so a subclass only has to override what it actually uses instead of the full interface.

On top of that, it adds:
- Its own `New()` runtime-selection (mesh, or mesh+dict+phase), mirroring OpenFOAM's `basicThermo::New<T>()` pattern.
- `Sct()` / `Prt()` (turbulent Schmidt/Prandtl numbers), set in the constructor from `constant/turbulenceProperties`: `Sct 0.4`/`Prt Sct` for LES, `0.7`/`Sct` for RAS, both `1.0` for laminar (`simulationType` missing defaults to laminar).

Calling an un-overridden virtual (e.g. `he()` on a class that never needed enthalpy) aborts via `NotImplemented`, not a silent no-op — a missing override shows up immediately at the call site.
