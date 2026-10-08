# Cantera Thermo
↑ [OxySim-129](./../../../README.md)  [basicThermo](./../basicThermoSTFS/README.md)

[canteraThermo.H](./canteraThermo.H)
[canteraThermo.C](./canteraThermo.C)
[canteraSpecieMixture.H](./canteraSpecieMixture.H)
[canteraSpecieMixture.C](./canteraSpecieMixture.C)

## Overview

`canteraThermo` implements `fluidThermo` (via [`basicThermoSTFS`](./../basicThermoSTFS)) and `basicSpecieMixture` (via `canteraSpecieMixture`) backed by a real Cantera `Solution`/`ThermoPhase`/`Transport` built from a mechanism file — mixture properties (`rho`, `psi`, `cp`, `cv`, `W`, `he`) and transport properties (`mu`, `nu`, `lambda`, per-species `D`/`Dtherm`) all come from Cantera rather than a lookup table.

```
mechanismFile "ch4_smooke_EKT.yaml";
phaseName     "gas";
energyType    ha;            // ha (absolute enthalpy) or ea (absolute energy)
transportModel UnityLewis;   // UnityLewis | Mix | Multi
useSoret      false;
inertSpecies  N2;
liquids { }
solids  { C; }
```

## Two roles in this repo

- **Standalone**: constructed with the default `initialUpdate = true`, so `updateTPY()`/`updateHPY()`/... actively drive all mixture/transport fields from Cantera — a normal, self-updating `fluidThermo`. Not currently used this way anywhere in this repo (both solvers use the second role below), but it's the default if you construct one directly.
- **Inside `flameletCloudFoam`/`flameletFoam`**: constructed with `initialUpdate = false`. None of its update functions are ever called in the solver's time loop — [`flameletThermo`](./../../flamelet/flameletThermo) owns `T`/`alpha`/`rho` instead (only `alpha` gets a manual one-way copy each PIMPLE iteration). Here `canteraThermo` only serves as a Cantera-backed species database (`speciesIndex`, `Wi(i)`, `Hc(i)`, element/atom counts via `canteraSpecieMixture`) for the particle submodels' stoichiometry.

## Testing

No dedicated unit test; exercised indirectly wherever an integration test constructs a `canteraThermo` alongside `flameletThermo`.
