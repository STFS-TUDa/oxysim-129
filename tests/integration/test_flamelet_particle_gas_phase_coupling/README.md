# Flamelet Particle-Gas Phase Coupling

**Author:** Pascal Steffens, **Date:** 29.01.2025

The goal of this test is to validate the flamelet model for Lagrangian cases. The test case consists of a particle in a channel that releases mass at a constant rate. This should lead to the correct mass at the outlet and the correct mass fraction of the released fuel. The test are carried out with up to three fuels.

## Description

The domain is a channel with dimensions 1 mm x 1 mm x 10 mm. The mesh consists of 3 x 3 x 100 cells. The gas flow is in z-direction. Boundary conditions are fixed value on the inlet, slip walls for velocity, and zero-gradient for all other quantitites.

![mesh](./reference/figures/mesh_particle_gas_phase_coupling.png)

## Cases

Three cases are tested with different fuel combinations.

### Case 0

Air at inlet, particle releases CH3, CH3O, and CH2O. 
Fuel mass fraction is supposed to reach 10 \% at outlet.

### Case 1

Air-methane mixture at inlet (5 \% methane). 
The particle releases CH3, CH3O, and CH2O. 
At the end of the domain, the released volatile mixture should reach 10 \%. This should lead to a reduction of the methane mass fraction to 4.5 \%. 

### Case 2 

Air-methane mixture at inlet (5 \% methane). 
The particle releases CH3, CH3O, and CH2O and CO. 
Oxygen is consumed from the gas phase to form CO.
At the end of the domain, the released fuel mixture should reach 8 \%, methane 4.5226 \%, and CO 3.52 \%.

## Flamelet table
The flamelet table is non-reactive with four dimensions: fuel mass fraction Z, blending factor 1 Y, blending factor 2 YY, and normalized enthalpy. The temperature range is 300 - 500 K. The fuels are a mixture of CH4, CH3O, and CH2O to represent volatiles, methane as a premixed assisting fuel, and char-off gases consisting of CO and N2. 

## Results
The figure shows the profiles for the transported quantities *fuel1*, *fuel2*, *fuel3*, and *ha*, and selected gas phase species and temperature. 
![results](reference/figures/plots.png)

## Notes
The tests are based on the theoretical values that the transported scalars and gas phase species should reach at the end of the domain. The calculations are in the `generate_reference_data.ipynb` notebook. A relative error of 0.5 \% is tolerated for all but CH3O. CH3O has a tolerance of 2 \%, as it has very small values in the case2 and therefore a very small absolute error would lead to failure of the test.

The present case can and should be extended to `reactiveCloudFoam`.
