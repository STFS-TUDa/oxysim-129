# Flamelet Equations
↑ [OxySim-129](./../../../README.md)

[EquationHandler.H](./EquationHandler.H)
[EquationHandler.C](./EquationHandler.C)
[inputVariableEquations.H](./inputVariableEquations.H)

## Overview

This directory holds the scalar transport equations for the flamelet input variables (e.g. mixture fraction `Z`, progress variable `Yc`, enthalpy `ha`) and the `EquationHandler` that registers and solves them.

Each variable listed under `variablesToSolve` in `constant/flameletProperties` is assigned an `updateType` in the `variableEquations` subdict, e.g.:

```
variableEquations
{
    Z
    {
        updateType unityLewis;
    }
}
```

`EquationHandler::appendFunctionMap()` maps each `updateType` string to a pair of functions: one to initialize the necessary fields, and one to build/solve the transport equation each time step. At runtime, `initializeAndRegisterEquations()` reads `variablesToSolve` and `variableEquations` from `flameletProperties` and registers the matching functions for every listed variable; `solve()` then invokes them in order.

Both `flameletFoam` and `flameletCloudFoam` use the same `EquationHandler` — pure-flamelet cases and particle-coupled cases just select different `updateType`s.

## Available `updateType`s

**Pure flamelet:**
- `noUpdate` — field is not solved (kept fixed)
- `unityLewis` — standard unity-Lewis-number scalar transport with a source term
- `unityLewis_noSource` — unity-Lewis-number transport without a source term
- `varianceAlgebraic` — algebraic model for a mixture fraction variance

**Particle–gas coupling (flameletCloudFoam):**
- `unityLewis_particleSource` — unity-Lewis-number transport with a particle-cloud source term
- `unityLewis_sourceTerm_particleSource` — as above, using a named source-term field
- `particleSource` — accumulates a particle-cloud source term into a field without solving a transport equation
- `unityLewis_particleEnthalpySource` — unity-Lewis-number enthalpy transport with a particle-cloud enthalpy source term
- `sum` — sets a field to the sum of other registered fields
- `fraction` — sets a field to the ratio of two other registered fields
- `scaleSourceTerm` — rescales a previously accumulated source term

## Adding a new equation

Each `updateType` is a pair of free functions in the `InputVariableEquations` namespace, with fixed signatures (`typedef`'d in `flameletTypeDefs.H`):

```cpp
// Called once at startup, to create any fields the update function needs
void InitFunctionType(std::string name, dictionary settings, flameletThermo& thermo);

// Called once per time step, to build and solve the transport equation (or do whatever else the updateType does)
void FunctionType(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions);
```

`name` is the variable's entry under `variablesToSolve` (e.g. `Z`); `settings` holds any extra keys from that variable's `variableEquations` subdict beyond `updateType` (empty if there are none). Fields aren't declared as local variables — the init function creates them by name on `thermo` (`thermo.appendField()` for a regular field, or `thermo.appendOutputField()` if it should additionally be looked up from the FLUT each step, e.g. a source term `omega_<name>`), and the update function retrieves them back by that same name via `thermo.getField()`.

To add a new one:
1. Declare the init/update function pair in `inputVariableEquations.H`.
2. Implement them in `inputVariableEquations.C` — `unityLewis`/`init_sourceTerm` is a good template for a standard transport equation with a FLUT source term; `noUpdate`/`init_empty` for a no-op.
3. Register the pair under a new `updateType` string in `EquationHandler::appendFunctionMap()`, e.g. `functionMap["myUpdateType"] = {init_myUpdateType, myUpdateType};`.
4. Add it to the "Available `updateType`s" list below.

`initializeEquations()` calls every registered variable's init function once (in `variablesToSolve` order) before `registerEquations()` wraps each update function into a closure; `solve()` then invokes those closures once per time step, in the same order. Both are keyed purely by the `updateType` string looked up in `functionMap` — an unregistered `updateType` fails at runtime (`std::invalid_argument`), not at compile time.
