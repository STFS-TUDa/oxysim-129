# Carbonaceous Cloud
↑ [OxySim-129](./../../../../README.md)

## Overview

Cloud-side counterpart to [`CarbonaceousParcel`](./../parcels/Templates/CarbonaceousParcel) — same layering as described in the top-level README's "OpenFOAM template architecture":

- `baseClasses/carbonaceousCloud/` — empty marker interface, only used for `isA<carbonaceousCloud>()` checks.
- `Templates/CarbonaceousCloud/` — adds devolatilization, surface-reaction and heat-capacity submodel handling on top of `CloudType`.
- `Templates/ThermoCloudSTFS/`, `ReactingCloudSTFS/`, `ReactingMultiphaseCloudSTFS/` — STFS's own forks of OpenFOAM's stock `ThermoCloud`/`ReactingCloud`/`ReactingMultiphaseCloud`, needed because the layers below `CarbonaceousCloud` must reference the STFS parcel submodel interfaces (`HeatCapacityModelSTFS`, `DevolatilisationModelSTFS`) instead of the stock ones.
- `derived/basicCarbonaceousCloud/` — the concrete typedef chain actually used by the solvers:

```cpp
CarbonaceousCloud<
    ReactingMultiphaseCloudSTFS<
        ReactingCloudSTFS<
            ThermoCloudSTFS<
                KinematicCloud<
                    Cloud<basicCarbonaceousParcel>>>>>>
```

`derived/basicReactingCloudSTFS/`, `basicReactingMultiphaseCloudSTFS/` and `basicThermoCloudSTFS/` are the same chain truncated at each intermediate level, instantiated for the corresponding STFS parcel type one level down.

See [`CarbonaceousParcel`](./../parcels/Templates/CarbonaceousParcel/README.md) for the parcel side and [`submodels`](./../submodels/README.md) for the devolatilization/heat-capacity/surface-reaction models plugged into this chain.
