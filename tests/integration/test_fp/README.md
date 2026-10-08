\ingroup GrpDummyPages
@defgroup flameletTesting
@ingroup flamelet

# Freely-propagating flame testcase 
**Author(s):** Matthias Steinhausen, Max Schneider, Pascal Steffens, Raphael Strickling
**Date:** 11.06.2024

## Description
In this testcase a 1D freely-propagating methane-air flame ($Z=0.05$) is calculated using tabulated chemistry.

## Initial and boundary conditions
The boundary conditions are configured as follows
| quantity | inlet (left) | outlet (right) |
|-|-|-|
| Velocity (U) | fixedValue (sl of table) | zeroGradient |
| Pressure (p) | zeroGradient | fixedValue |
| Mixture Fraction (Z) | fixedValue (freshGas) | zeroGradient |
| Progress Variable (PV) | fixedValue (freshGas) | zeroGradient |


The progress variable and the velocity is initialized from a freely propagating Cantera flame. For this, the flame is sampled at 10 points and then mapped to the $x$-coordinate in OpenFOAM. In this case, the progress variable is set to $PV = CO + CO_\text{2}$

![Initial profiles for the progress variable and the velocity. The velocity at the inlet on the left side is the laminar flame speed.](reference/figures/init_profiles.png){width=75%}


The mixture fraction $Z$ is initialized with a uniform value of $Z = 0.05$

## Quantities of interest
In the testcase different quantities are compared to a reference OpenFOAM solution:
+ Global flame properties (Laminar flame thickness, displacement speed, consumption speed, flame position)
+ OpenFOAM solution variables (PV, Z, U, p, omega_yc and rho)
+ Maximum and Minimum value
    + Integral under the curve
    + Area between the curves (normalized by the integral of the curve)
    + Pearson correlation coefficient
+ Size of the pressure fluctuations

## Results
Results from the last test with the `--write_report` True flag:

<figure>
    <img src="report/figures/t_0.01_physicalSpace.png" alt="Figure is missing. Did you run the testing with the write_report flag?" width="800"/>
    <figcaption>Comparison in physical space with OpenFOAM reference.</figcaption>
</figure>
<figure>
    <img src="report/figures/t_0.01_stateSpace.png" alt="Figure is missing. Did you run the testing with the write_report flag?" width="800"/>
    <figcaption>Comparison in state space with OpenFOAM reference.</figcaption>
</figure>
<details>
    <summary>Reference pictures how it should look like</summary>
    ![Comparison in physical space with OpenFOAM reference.](reference/figures/t_0.01_physicalSpace.png){width=75%}
    ![Comparison in state space with OpenFOAM reference.](reference/figures/t_0.01_stateSpace.png){width=75%}
</details>
<details>
  <summary>Comparison of the reference with cantera</summary>
    ![Comparison of the OpenFOAM reference solution with cantera.](reference/figures/compare_openfoam_cantera.png){width=75%}
</details>