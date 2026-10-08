# Flamelet Thermo
↑ [OxySim-129](./../../../README.md)    [basicThermo](./../../thermo/basicThermoSTFS/README.md)

[flameletThermo.H](./flameletThermo.H)
[flameletThermo.C](./flameletThermo.C)

## Overview

`flameletThermo` implements OpenFOAM's `fluidThermo` interface (via [`basicThermoSTFS`](./../../thermo/basicThermoSTFS)) by looking up thermochemical quantities from a tabulated flamelet look-up table (FLUT), read with [`flut-reader`](https://github.com/STFS-TUDa/flut-reader)'s `FLUT::LookupTable` (`extern/flut-reader`). Both `flameletFoam` and `flameletCloudFoam` use it as their thermo object.

On `update()` it reads the current cell/face **input variables** (e.g. `Z`, `yc`, `ha`), looks them up in the FLUT, and writes the resulting **output variables** (T, density, viscosity, species mass fractions, ...) into shared `volScalarField`s.

## `flameletProperties` dictionary

```
flutFile                     FLUT.h5;              // FLUT file, default: FLUT.h5
inputVariablesFlameletTable  (Z yc);                // must match the FLUT's own input variables exactly
outputVariablesFlameletTable (N2 O2 CH4);           // additional fields to look up, beyond the built-in defaults
variablesToSolve              (Z yc);               // order of the transport-equation variables (see flameletEquations)
useWithCanteraThermo         false;                 // see note below — leave at false
```

- `inputVariablesFlameletTable` is just a consistency check against the FLUT's own input list; a mismatch is a `FatalError`.
- With the default `useWithCanteraThermo false`, the default output variables (see [Required FLUT outputs](#required-flut-outputs) below) are refreshed from the FLUT every `update()`. `thermo:psi` is computed manually afterwards as `thermo:rho / p` (or `/ pRef`), not looked up. Anything extra goes in `outputVariablesFlameletTable`.
- **`useWithCanteraThermo true` is effectively dead**: it's meant to let `canteraThermo` own `T`/`alpha`/`mu` instead, but that branch never registers them anywhere, so they just keep their old values (`mu()` even throws). No solver sets this flag, and `flameletCloudFoam`'s `canteraThermo` is never `correct()`ed in the time loop anyway — only the unit test exercises this path, and it asserts fields stay untouched. Don't rely on it.
- A field can't be both a `variablesToSolve` entry and an output/input variable — enforced with a `FatalError`.
- [`flameletEquations`](./../flameletEquations) solves `variablesToSolve`; `flameletThermo` handles everything else via table lookup.

## Required FLUT outputs

What the FLUT must actually supply as output variables (i.e. its `output_dict` entries) depends on which solver and models are in play — this pulls together requirements that otherwise only show up by reading the source:

| Output | Required when | Where it's consumed |
|---|---|---|
| `T`, `alpha`, `rho`, `mu` | Always, for any flamelet case | `flameletThermo::init_outputVariables()` (literal names, not configurable) |
| `cp` | `flameletCloudFoam` only | appended in `createFields.H`, used by `Cpv()` |
| Every species in the Cantera mechanism that isn't itself a `variablesToSolve` entry | `flameletCloudFoam` only | `createFields.H` loops `thermo.species()` and appends whichever aren't already transported |
| `W` (mean molecular mass) | Radiation model active (`wsggBordbarAbsorptionEmission`/`wsggDorigonAbsorptionEmission`) | `flameletThermo::W()`, called from the WSGG absorption-emission models — see [`canteraThermo`](./../../thermo/canteraThermo) and [`radiation`](./../../radiation) |
| `ha` in the **boundary** FLUT | A `boundaryLookUp` BC is used | see [Non-adiabatic boundary lookup](#non-adiabatic-boundary-lookup) below |
| `yc_max` in the **bulk** FLUT | `specialBoundaryLookup` only | same |

A missing output fails loudly rather than silently: `flut-reader`'s `OutputDict::containsVariable` throws `Requested outputVariable: '<name>' not present in FLUT.` the first time it's looked up — typically at the first `update()` call, not at construction. Anything beyond this list (e.g. species kept only for post-processing) is a free choice via `outputVariablesFlameletTable`.

`flut-reader`'s [demo notebook](https://github.com/STFS-TUDa/flut-reader/blob/main/flut_creation_demo/flamelet_demo.ipynb) builds a table from Cantera flame solutions and annotates each exported column against this same list.

## Non-adiabatic boundary lookup

An `ha` boundary of type `boundaryLookUp` ([`flameletBCs/boundaryLookUp`](./../flameletBCs)) triggers a second, boundary-only FLUT lookup from `flutFileBC` (default `FLUT_BC.h5`), writing `ha` only at the flagged patches. All such patches in a case must use the same `lookupType`:
- `standardBoundaryLookup` — boundary FLUT queried directly with the patch's own inputs (default).
- `specialBoundaryLookup` — legacy Ketelheun (2013) treatment, adds a normalized `ccBoundary = yc / yc_max` input.

## Field ownership (`appendField` / `getField`)

Output fields aren't fixed members — `flameletThermo` keeps a `std::map<word, volScalarField&>` (`fields_`) built on the fly:

- `appendField(...)` reuses the field from the object registry if it already exists (e.g. created by `canteraThermo` or the solver), otherwise constructs and owns a new one.
- `appendOutputField(...)` = `appendField(...)` + registers it for FLUT lookup on every `update()`.
- `getField(name)` is the only accessor; missing names raise a `FatalError`.

This is what lets other components (e.g. the radiation model's `W()` call) share fields with `flameletThermo` without a compile-time dependency on who created them.

## Testing

Unit-tested by `tests/unit/testFlameletThermo` (Catch2, `COMPILE_LAGRANGIAN_WITH_TESTS=1`) — see the top-level README's "Testing" section for how to build/run it.
