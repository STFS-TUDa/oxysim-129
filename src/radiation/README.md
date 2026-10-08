# Radiation
↑ [OxySim-129](../../README.md)

## Overview

This library adds a weighted-sum-of-gray-gases (WSGG) radiation model on top of OpenFOAM's finite-volume discrete ordinates method (fvDOM), ported from `reactivefoam_lagrangian_radiation`. It is used by `flameletCloudFoam` to solve the radiative transfer equation (RTE) for a non-grey participating medium and to couple the resulting source term into the gas- and particle-enthalpy equations.

```
radiation/
├── radiationModel/
│   ├── radiationModelSTFS/       Top-level model: resolves the basicThermo to radiate against
│   └── fvDOMwsgg/
│       ├── fvDOMwsgg/            Finite-volume discrete-ordinates solver (non-grey, WSGG)
│       └── radiativeIntensityRaywsgg/  Single discrete-ordinate ray intensity field
├── submodels/
│   ├── wsggAbsorptionEmission/           Base class: gas-phase WSGG absorption/emission
│   ├── wsggBordbarAbsorptionEmission/    WSGG coefficients fitted for oxyfuel atmospheres
│   ├── wsggDorigonAbsorptionEmission/    WSGG coefficients fitted for air atmospheres
│   └── wsggBinaryAbsorptionEmission/     Combines a gas-phase WSGG model with a second
│                                         absorption model (e.g. a particle cloud)
└── derivedFvPatchFields/
    └── wideBandDiffusiveRadiationwsgg/   Wall BC: wide-band diffusive radiation at a given
                                          patch temperature
```

## `fvDOMwsgg`

`fvDOMwsgg` solves the RTE for `nPhi × nTheta` discrete directions per octant, without scatter or reflective walls, using a non-grey WSGG representation of the gas mixture (see Trivic, doi:[10.1016/j.ijheatmasstransfer.2003.09.027](https://doi.org/10.1016/j.ijheatmasstransfer.2003.09.027) for the underlying formulation). It is selected as `radiationModel` in `constant/radiationProperties`:

```
radiation       on;
radiationModel  fvDOMwsgg;

fvDOMwsggCoeffs
{
    nPhi        2;      // azimuthal angles in PI/2 on X-Y (from Y to X)
    nTheta      4;      // polar angles in PI (from Z to X-Y plane)
    maxIter     5;      // maximum number of iterations per radiation solve
    tolerance   1e-3;   // convergence tolerance for the radiation iteration
}

solverFreq 1;            // flow iterations per radiation iteration
thermoName flameletProperties;  // basicThermo object to radiate against
```

`thermoName` disambiguates which `basicThermo` object `radiationModelSTFS` should look up when more than one is registered (e.g. `flameletProperties` for the flamelet gas-phase thermo).

## Gas-phase absorption/emission (WSGG)

`wsggAbsorptionEmission` is the common base for gray-gas weighting-factor models; it is not used directly but selected through one of its two coefficient sets:

- **`wsggBordbarAbsorptionEmission`** — fitted for oxyfuel (CO₂/H₂O-rich) atmospheres. Bordbar et al. (2014), doi:[10.1016/j.combustflame.2014.03.013](https://doi.org/10.1016/j.combustflame.2014.03.013)
- **`wsggDorigonAbsorptionEmission`** — fitted for air atmospheres. Dorigon et al. (2013), doi:[10.1016/j.ijheatmasstransfer.2013.05.010](https://doi.org/10.1016/j.ijheatmasstransfer.2013.05.010)

Each carries its own polynomial coefficient set (`c1Coeffs`…`dCoeffs` for Bordbar, `biCoeffs`/`kiCoeffs` for Dorigon) in the corresponding `<model>Coeffs` subdict — see the header of each class for the full dictionary example.

## Coupling particle radiation into the gas phase

To combine the gas-phase WSGG model with radiation from a Lagrangian carbonaceous cloud, select `wsggBinaryAbsorptionEmission` as the top-level `absorptionEmissionModel` and nest the two contributions under `model1`/`model2`:

```
absorptionEmissionModel wsggBinaryAbsorptionEmission;
wsggBinaryAbsorptionEmissionCoeffs
{
    model1
    {
        absorptionEmissionModel wsggDorigonAbsorptionEmission;
        wsggDorigonAbsorptionEmissionCoeffs { ... }
    }
    model2
    {
        absorptionEmissionModel wsggGreyCloudAbsorptionEmission;
        wsggGreyCloudAbsorptionEmissionCoeffs
        {
            cloudNames ( coalCloud ashCloud );
        }
    }
}
```

`model1` must be a gas-phase WSGG model (Bordbar or Dorigon); `model2` provides the particle-cloud contribution. The particle side of this coupling — `wsggGreyCloudAbsorptionEmission` — lives with the other Lagrangian submodels; see [Submodels → Cloud radiation](../lagrangian/intermediate/submodels/README.md#cloud-radiation).

## Wall boundary condition

`wideBandDiffusiveRadiationwsgg` provides a wide-band diffusive radiation condition for a patch at a given temperature, with one emissivity per WSGG band:

```
".*"
{
    type            opaqueDiffusive;
    wallAbsorptionEmissionModel
    {
        type multiBandAbsorption;
        emissivity      (1 1 1 1 1);   // one value per WSGG band (typically 4 grey + 1 clear)
    }
}
```

## Solver integration

`flameletCloudFoam` calls `radiation->correct()` once per PIMPLE iteration (`flameletCloudFoam.C`), right after the flamelet table update. The resulting energy source term is added directly into the gas-enthalpy transport equation through the `unityLewis_particleEnthalpySource` `updateType` in [`EquationHandler`](../flamelet/flameletEquations/README.md) (`inputVariableEquations.C`), which sums the particle-cloud enthalpy source and `radiation.Sh(thermo, controlVariable)` into the same matrix. Particle-side radiative exchange (`radiation semiImplicit 1;` under a cloud's `sourceTerms` dict, `radiation on;` in the cloud properties) is handled by OpenFOAM's Lagrangian radiation coupling together with `wsggGreyCloudAbsorptionEmission`.

Required libraries at run time (see `system/controlDict` → `libs`): `libradiationModelsSTFS.so` and `libparticleRadiation.so`.

## Testing

- [`test_multiphaseRadiation_dorigon`](../../tests/integration/test_multiphaseRadiation_dorigon/README.md) — air-atmosphere case, `wsggDorigonAbsorptionEmission` + `wsggGreyCloudAbsorptionEmission`
- [`test_multiphaseRadiation_bordbar`](../../tests/integration/test_multiphaseRadiation_bordbar/README.md) — oxyfuel-atmosphere case, `wsggBordbarAbsorptionEmission` + `wsggGreyCloudAbsorptionEmission`, with additional particle-emissivity, area-density and temperature-offset sweeps
