# Multiphase Radiation Test — Dorigon (air atmosphere)

This test validates the coupled gas + particle radiation model (see [Radiation](../../../src/radiation/README.md) and [Cloud radiation](../../../src/lagrangian/intermediate/submodels/README.md#cloud-radiation)) for an air atmosphere against the reference case of Gronarz (2017).

## Description

A fixed cloud of coal and ash particles (pre-populated Lagrangian restart data in `0.org/lagrangian/{coalCloud,ashCloud}/`, injection disabled in both `coalCloudProperties`/`ashCloudProperties`) is placed in a radiating, non-scattering gas medium. `fvDOMwsgg` solves the RTE with `wsggDorigonAbsorptionEmission` (air-atmosphere WSGG coefficients) for the gas phase, combined via `wsggBinaryAbsorptionEmission` with `wsggGreyCloudAbsorptionEmission` for the particle-cloud contribution.

## Test

Three cases sweep the particle emissivity `epsilon0` in `constant/coalCloudProperties`:

| Test | `epsilon0` |
|---|---|
| `test_eps0` | 0 (no particle emission/absorption) |
| `test_eps029` | 0.29 (default) |
| `test_eps1` | 1 (black particles) |

For each case, the radiative source term `Sh` (gas absorption/emission plus the particle emission/absorption terms computed from the parcels' `radAreaP`/`radAreaPT4` fields) is sampled along a line through the domain and numerically compared (`np.testing.assert_allclose`, `rtol=1e-2`, `atol=0.01`) against `reference/foam/epsi_coal*.parq`, a previously verified OpenFOAM result (regression reference). `reference/Gronarz/epsi_coal*.parq` — the literature reference case of Gronarz (2017) — is plotted alongside it in the optional `--write_report` figure, but is not part of the pass/fail assertion.

## Reference solution

![epsilon sweep](./reference/figures/epsi_coal.png)
