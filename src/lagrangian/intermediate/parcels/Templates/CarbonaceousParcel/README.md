# Carbonaceous Cloud and Parcel
↑ [OxySim-129](../../../../../../README.md)

[CarbonaceousParcel.H](./CarbonaceousParcel.H)  
[CarbonaceousCarcel.C](./CarbonaceousCarcel.C)  
[CarbonaceousCarcel.C](./CarbonaceousParcelI.H)  
[CarbonaceousCarcelIO.C](./CarbonaceousParcelIO.C)  

## Basics

The carbonaceous cloud and parcel are used to simulate coal and biomass particles. Submodels for carbonaceous parcels and other parcels are described [here](../../../submodels/README.md).

The CarbonaceousParcel template is derived from the ReactingMultiphaseParcel. The inheritance structure is visible in [basicCarbonaceousParcel.H](../../derived/basicCarbonaceousParcel/basicCarbonaceousParcel.H). 

Notable differences are the inclusion of the DevolatilisationModelSTFS and the HeatCapacityModelSTFS. The DevolatilizationModelSTFS is a special version of the OpenFOAM DevolatilisationModel class that handles also the mass that is converted to char. The HeatCapacitySTFS model does not have a counterpart in the original OF. The SurfaceReactionModelSTFS has the same inputs as the original OF version, whose update function is called in the original OF MultiphaseReactingParcel.

## Notable differences to the MultiphaseReactingParcel

Compared to the MultiphaseReactingParcel, the CarbonaceousParcel has the *dMass_ipc_fromGas* and *dMass_ipc_toSolid* variables, with ipc for internal phase change. For every phase change (only gas to solid currently implemented) the *to* and *from* variable is needed. The devolatilization model makes use of this to convert some of the gas mass to solid char. The mass exchange is calculated in the devolatilization model. But only in the *updateMassFractions* function the terms are added in the arguments. The function expects positive numbers to be a loss. Therefore, the fromGas variable is added and toSolid is subtracted. 