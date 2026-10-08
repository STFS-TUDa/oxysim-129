# Tabulated Chemistry Particle Heating Test

**Author:** Pascal Steffens **Date:** 22.01.2025

The test aims to verify that the particle sees the correct gas phase temperature and heats up to the same temperature, using tabulated chemistry. 
The interpolation method for temperature that is chosen here is cellPoint as this has shown to be error prone.

*In the future the test could be updated to include testing of the heat transfer.*

## Description
The test case is carried out in a small cubical domain representing a channel section with dimensions 1 mm x 1 mm x 1 mm that is shown in the figure. 
The mesh consists of eight cells.
The **wall boundary** has a **slip** condition for velocity and **zero-gradient** for the remaining variables. 
A pure nitrogen stream enters through the inlet in positive z-direction with a velocity of 0.4 m/s and exits through the outlet patch. 

![mesh](./reference/figures/mesh_heat_transfer_test.png)

### Particles

Three particles are positioned in one cell, at a corner point, between the corner and the cell center, and at the cell center (see mesh). 
Particle | x (mm) | y (mm) | z (mm) |
| - | - | - | - | 
| 1 | 0.5 | 0.5 | 0.5 |
| 2 | 0.625 | 0.625 | 0.625 |
| 3 | 0.75 | 0.75 | 0.75 |

Particles are set up to couple only in one direction, meaning that the gas phase is unaffected by the particle heating. 
This allows validating that the particle correctly equilibrates to the gas phase temperature.

### Temperature profiles
The test case runs for 0.07 s. 
The temperature is adapted by applying a **uniformFixedValue** boundary condition for enthalpy to the inlet that allows to define a time dependent profile.

| Time (s) | Temperature (K) | Enthalpy (J/kg) |
| - | - | - |
| 0 - 0.01 | 300 | 1971 |
| 0.01 - 0.04 | 425 |133000 |
| 0.04 - 0.07 | 675 |400000 | 

### Flamelet
The flamelet table has a single dimension along absolute enthalpy `ha` (`inputVariablesFlameletTable (ha)`) and only contains nitrogen. Therefore, only the enthalpy equation is solved.

## Reference solution
The reference solution shows the temperature over time for probes and the particle cloud. 
Since the enthalpy change enforced on the inlet, the gas phase temperature profiles have no sharp edges. 
Further, the temperature change first reaches the first position, which is at the smallest z-coordinate, and the third position last.
Particle temperatures follow. 
It is important that particles eventually reach the exact same temperature as the gas phase.

![reference_profiles](./reference/figures/reference_profiles.png)


 
