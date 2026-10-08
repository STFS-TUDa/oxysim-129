# Submodels
↑ [OxySim-129](../../../../README.md)  
This is an overview of the available custom submodules.

## Carbonaceous
### Devolatilization Model STFS
#### [C2SM](./Carbonaceous/DevolatilizationModelSTFS/C2SM/C2SM.H)
Kobayashi Competing Two Step Model
It describes the reaction of raw coal or biomass to char and volatiles that are released to the environment. The model was presented by [Kobayashi et al. (1977)](https://doi.org/10.1016/S0082-0784(77)80341-X). The implementation is visibile in the [code](./Carbonaceous/DevolatilizationModelSTFS/C2SM/C2SM.C) and usage is shown in the [header](./Carbonaceous/DevolatilizationModelSTFS/C2SM/C2SM.H). In [Steffens et al. (2024)](https://doi.org/10.1016/j.fuel.2024.132098) it is introduced as follows.
> The devolatilization of coal and biomass is calculated with a competing two-step model based on Kobayashi [[35](https://doi.org/10.1016/S0082-0784(77)80341-X)]. The conversion of particle mass to volatiles and char is defined with the following equation 
> $$ Y(t) = \int_0^t (Y_1 k_1 + Y_2 k_2) \exp \left( - \int_0^{t'} (k_1 + k_2) \mathrm{d} \tau \right) \mathrm{d} t' $$
> with $Y_1$ and $Y_2$ corresponding to the final volatile yield of the low and high heating rate reactions, respectively, and the equations $k_1$ and $k_2$
> $$k_i(T_p) = A_i T_P \exp\left(-\frac{E_{A,i}}{R_G T_P}\right)$$

The parameters of the equations can be obtained from detailed solid chemistry simulations using the Pyrolysis Kinetics Pre-Processor ([pkp-lite](https://github.com/STFS-TUDa/pkp-lite)). This approach was described by [Vascellari et al. (2013)](http://dx.doi.org/10.1016/j.fuel.2013.06.014).

The C2SM model is [tested](../../../../tests/integration/test_C2SM/test_C2SM.py), the [notebook](../../../../tests/integration/test_C2SM/generate_reference_data.ipynb) for reference data gives an idea of the behavior.

### Heat Capacity Model STFS
#### [ConstantHeatCapacity](./Carbonaceous/HeatCapacityModelSTFS/ConstantHeatCapacity/ConstantHeatCapacity.H)
The particle keeps a constant heat capacity.
#### [HeatCapacityFit](./Carbonaceous/HeatCapacityModelSTFS/HeatCapacityFit/HeatCapacityFit.H)
The heat capacity changes according to the temperature with a polynomial function.
$$ c_p = a + b T + c T^2 + d T^3 + e T^4$$
From [Steffens et al. (2024)](https://doi.org/10.1016/j.fuel.2024.132098)
> The heat capacity is obtained using a polynomial of the particle temperature [[33](http://dx.doi.org/10.1021/acs.energyfuels.0c03479),[34](http://dx.doi.org/10.1016/0016-2361(96)00067-1)]. Given the limited availability of specific heat capacity data for walnut shells, the heat capacity polynomial of lignite is adopted for both fuels. Both coal and biomass are pre-dried.  
### Surface Reaction Model STFS
#### [COxidationKineticDiffusionLimitedRateSTFS](./Carbonaceous/SurfaceReactionModelSTFS/COxidationKineticDiffusionLimitedRateSTFS/COxidationKineticDiffusionLimitedRateSTFS.H)

This model is basically the Baum-Street-Model, which describes the oxidation of char. In the present version, the char is oxidized to CO. The original source is [Baum and Street (1971)](https://doi.org/10.1080/00102207108952290).
In [Steffens et al. (2024)](https://doi.org/10.1016/j.fuel.2024.132098) it is described as follows:
>Char conversion is calculated with the model by Baum and Street [[39](https://doi.org/10.1080/00102207108952290)], which uses the following correlation to describe the conversion of the char structure
>$$ \frac{\mathrm{d} m_{\mathrm{char}}}{\mathrm{d} t}=-\pi D_{\mathrm{prt}}^2 \rho R T \frac{Y_{\mathrm{O_2}}}{W_{\mathrm{O_2}}}\left(R_{\mathrm{diff}}^{-1}+R_{\mathrm{reac}}^{-1}\right)^{-1} $$
>with $D_{\mathrm{prt}}$ as particle diameter, $Y_{\mathrm{O_2}}$ oxygen mass fraction, and $W_{\mathrm{O_2}}$ oxygen molar mass. The diffusion rate coefficient is defined as
>$$R_{\mathrm{diff}}=\frac{C_{\mathrm{diff}}}{D_{\mathrm{prt}}}\left(\frac{T_{\mathrm{prt}}+T}{2}\right)^{0.75}.$$
>The reaction rate coefficient $R_{\mathrm{reac}}$ amounts to 
>$$ R_{\mathrm{reac}}=A\exp\left(\frac{-E_a}{RT_{\mathrm{prt}}}\right).$$
>The coefficients $C_{\mathrm{diff}}$, $A$, and $E_a$ were determined experimentally by Gövert et al. [[40](http://dx.doi.org/10.1016/j.fuel.2017.03.009)] and are shown in Table 3.
>Biomass particles, being larger, have a diminished influence on the gas phase due to their significantly lower fixed carbon content compared to coal.
>Studies have shown that the combustion of lignocellulosic chars is similar to coal char combustion [[41](http://dx.doi.org/10.1016/j.pecs.2008.08.001)], leading to the assumption that both coal and biomass chars are essentially pure carbon [[7](http://dx.doi.org/10.1155/2018/7036425),[22](http://dx.doi.org/10.1016/j.proci.2022.07.135)].
>The heat of the heterogeneous reaction is given to the particle [[15](http://dx.doi.org/10.1016/j.combustflame.2016.07.013),[32](http://dx.doi.org/10.1016/j.fuel.2020.117683)].

## Cloud radiation
### [wsggGreyCloudAbsorptionEmission](./CloudRadiation/wsggGreyCloudAbsorptionEmission/wsggGreyCloudAbsorptionEmission.H)

Extends the gas-phase WSGG [radiation model](../../../radiation/README.md) with radiation from Lagrangian carbonaceous particle clouds. It computes the volumetric particle emissive power from parcel temperatures, surface areas and a constant particle emissivity, distributes it into the WSGG bands using gas weighting factors evaluated at the particle temperature, and forms a gray volumetric particle absorption coefficient from the time-averaged particle surface area per cell.

It is selected as `model2` of a [`wsggBinaryAbsorptionEmission`](../../../radiation/README.md#coupling-particle-radiation-into-the-gas-phase) model, alongside a gas-phase WSGG model (Bordbar or Dorigon) as `model1`:

```
absorptionEmissionModel wsggGreyCloudAbsorptionEmission;
wsggGreyCloudAbsorptionEmissionCoeffs
{
    cloudNames
    (
        coalCloud
    );
}
```

Tested by [test_multiphaseRadiation_dorigon](../../../../tests/integration/test_multiphaseRadiation_dorigon/README.md) and [test_multiphaseRadiation_bordbar](../../../../tests/integration/test_multiphaseRadiation_bordbar/README.md).

