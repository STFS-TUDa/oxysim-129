# Multiphase Radiation Test — Bordbar (oxyfuel atmosphere)

This test validates the coupled gas + particle radiation model (see [Radiation](../../../src/radiation/README.md) and [Cloud radiation](../../../src/lagrangian/intermediate/submodels/README.md#cloud-radiation)) for an oxyfuel atmosphere against the reference case of Gronarz (2017).

## Description

A fixed cloud of coal and ash particles (pre-populated Lagrangian restart data in `0.org/lagrangian/{coalCloud,ashCloud}/`, injection disabled in both `coalCloudProperties`/`ashCloudProperties`) is placed in a radiating, non-scattering gas medium. `fvDOMwsgg` solves the RTE with `wsggBordbarAbsorptionEmission` (oxyfuel-atmosphere WSGG coefficients) for the gas phase, combined via `wsggBinaryAbsorptionEmission` with `wsggGreyCloudAbsorptionEmission` for the particle-cloud contribution.

This case exercises a broader parameter sweep than [test_multiphaseRadiation_dorigon](../test_multiphaseRadiation_dorigon/README.md):

| Test | Parameter swept |
|---|---|
| `test_standard` | default settings |
| `test_eps0` / `test_eps1` | particle emissivity `epsilon0` = 0 / 1 |
| `test_Ap025` / `test_Ap4` | coal particle diameter `d` (1e-2 m / 4e-2 m vs. the 2e-2 m baseline), which changes the resulting projected-area concentration `Ap` |
| `test_deltaT0` / `test_deltaT200` | gas–particle temperature offset `deltaT` |

For each case, the radiative source term `Sh` (gas absorption/emission plus the particle emission/absorption terms computed from the parcels' `radAreaP`/`radAreaPT4` fields) is sampled along a line through the domain and numerically compared (`np.testing.assert_allclose`, `rtol=1e-2`, `atol=0.01`) against `reference/foam/*.parq`, a previously verified OpenFOAM result (regression reference). `reference/Gronarz/*.parq` — the literature reference case of Gronarz (2017) — is plotted alongside it in the optional `--write_report` figure, but is not part of the pass/fail assertion.

## Reference solution

![epsilon sweep](./reference/figures/epsi_coal.png)
![particle area density sweep](./reference/figures/Ap.png)
![temperature offset sweep](./reference/figures/deltaT.png)
