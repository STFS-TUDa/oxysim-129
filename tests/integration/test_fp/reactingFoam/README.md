# reactingFoam reference case

Same flame as the flameletFoam test case (`../foam_template`), computed with the standard OpenFOAM v2512 solver `reactingFoam` and detailed chemistry.
Not run during testing; the result is stored in `../reference/foam/reactingFoam_t_0.010.parq` and shown in the report figures and in `../reference/figures/compare_openfoam_cantera.png`.

## Setup
| | flameletFoam | reactingFoam |
|-|-|-|
| Chemistry | FLUT (Z, yc) | Smooke mechanism (16 species, 25 reactions), `chemkinToFoam`, operator splitting with seulex (absTol 1e-12, relTol 1e-6) |
| Species diffusion | unity Lewis | unity Schmidt (`muEff`) |
| Transport | from FLUT | Sutherland, air coefficients for all species (`chemkin/transportProperties`) |
| Pressure coupling | low Mach | compressible |
| Outlet pressure | fixedValue | waveTransmissive (non-reflecting; required for the compressible solver) |

Identical in both cases: mesh (0.025 m, 2500 cells), time step (1e-6 s, end time 0.01 s), PIMPLE (1 outer, 3 inner correctors), discretization schemes, initial conditions (`initialize.py`: species and T from `../reference/cantera/fp_Z0.05000.parq`, sampled as for flameletFoam; uniform U = $S_L$) and all other boundary conditions.
flameletFoam gives the same result with a waveTransmissive outlet.

## Run
```
module load openfoam/2512
./Allrun -p
```

## Results (t = 0.01 s)
Deviation from the Cantera reference `../reference/cantera/fp_Z0.05000.parq` (unity Lewis), computed in `../create_and_validate_reference_solutions.ipynb`. Compute time: solver ExecutionTime, 5 MPI processes on one exclusive node (2x Intel Xeon E5-2680 v2, 2.8 GHz).

| | flameletFoam | reactingFoam |
|-|-|-|
| Consumption speed $S_c$ | -2.8 % | +0.4 % |
| Displacement speed $S_d$ | -2.2 % | +2.0 % |
| Thermal thickness $\delta_T$ | -2.2 % | -1.8 % |
| Peak $\dot\omega_{yc}$ | -0.2 % | +24.5 % |
| Burnt gas temperature $T_b$ | 0 K | -11 K |
| Compute time | 34 s | 2038 s |

Global flame properties agree with the reference within 3 % for both solvers. In state space, flameletFoam reproduces the reference, whereas reactingFoam has a higher T and lower rho at a given yc and a higher $\dot\omega_{yc}$ peak shifted to lower yc. Potential causes: reactingFoam diffuses species with unity Schmidt number, whereas the reference and the FLUT use $Le = 1$.

The chemistry integration settings of reactingFoam (operator splitting, ODE tolerances) are not optimized for speed; independent of this, reactingFoam deviates from the reference in state space.

![Comparison of flameletFoam and reactingFoam with the Cantera reference in state space.](../reference/figures/compare_openfoam_cantera.png)
