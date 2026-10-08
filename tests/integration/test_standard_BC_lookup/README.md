# Standard Boundary-Lookup Test

Validates the `standardBoundaryLookup` path of [`flameletThermo`](./../../../src/flamelet/flameletThermo)'s non-adiabatic boundary treatment (as opposed to `specialBoundaryLookup`, covered by the `testFlameletThermo` unit test).

## Setup

Not a physically meaningful flow case: the `WALL` patch's 100 boundary faces are each initialized with an independent, effectively random point in the flamelet input space (`Z`, `yc`, `ha` in `0.org/{Z,yc,ha}`, spanning close to their full ranges). `ha` there is a `boundaryLookup` BC with `lookupType standardBoundaryLookup`, so `flameletThermo` looks each face up in the boundary FLUT (`FLUT_BC.h5`, via `flutFileBC`) and writes the resulting bulk-table temperature `T`. A second, independently-prescribed field `T_BC` (`0.org/T_BC`, `fixedValue` on `WALL`) holds the expected temperature for each of those 100 points, generated alongside `FLUT.h5`/`FLUT_BC.h5` in `create_flut.ipynb`.

`flameletFoam` runs for 100 steps (`deltaT 1e-5`, `endTime 1e-3`) so the boundary lookup settles once; `variablesToSolve` (`Z`, `ha`, `yc`) all use `updateType unityLewis_noSource`.

## Test

`test_standard_BC_lookup.py` reads `T` and `T_BC` on the `WALL` patch at `t = 0.001` (via `pyvista`) and asserts `T ≈ T_BC` (`atol=10` K) — i.e. that the boundary FLUT lookup reproduces the reference temperature at all 100 sampled points. With `--write_report`, the per-face values and error are written to `report/report.csv`.
