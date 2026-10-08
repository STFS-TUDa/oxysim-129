# Competing Two-Step Model Test

**Author:** Pascal Steffens **Date:** 12.02.2025

This test verfies the Competing-Two-Step-Model (C2SM model), which computes the conversion to volatiles and char. 
This model was introduced by [Kobayashi 1977](https://doi.org/10.1016/S0082-0784(77)80341-X) (find the formulation there). 
An implementation can be found in the C2SM devolatilization submodel in carbonaceous parcel submodels and in the `generate_reference_data.ipynb` notebook.

## Description

The domain is the same as in the particle-gas phase coupling test case (see [here](../test_flamelet_particle_gas_phase_coupling/README.md)).

### Boundary conditions
The temperature is set to 800 K. Inflow velocity of pure nitrogen is 1 m/s.

### Flamelet
The flamelet table has two dimensions along mixture fraction `Z` and absolute enthalpy `ha` (`inputVariablesFlameletTable (Z ha)`). Equations for particle source term, Z, and enthalpy are solved. 

### Particle

The domain contains a parcel composed of 10 particles at z = 5 mm. The diameter is 30 µm and the combined mass 1e-10 kg. At the end the particles will have shrunk to ca. 20 µm and lost 70 \% of their mass. Therefore, particles have a final mass of ca. 3e-12 kg and the volatiles that pass the outlet should have a mass of 7e-11 kg.

## Test

For the test, the C2SM results for mass release over time, char mass fraction over time, and volatiles in gas phase are calculated in the `generate_reference_data.ipynb` notebook. 

The test compares the time resolved data for char mass fraction and particle mass, and the final particle mass and total volatile yield at the end of the domain.

## Reference solution

The reference solution is calculated in python with an explicit solution algorithm. 

Results python:

![Python](./reference/figures/C2SM_python.png)

Results OpenFOAM:

![OpenFOAM](./reference/figures/C2SM_openfoam.png)

Comparison:

![comparison](./reference/figures/comparison.png)



 
