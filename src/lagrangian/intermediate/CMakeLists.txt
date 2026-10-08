set(CLOUDS    ${CMAKE_CURRENT_SOURCE_DIR}/clouds)
set(PARCELS   ${CMAKE_CURRENT_SOURCE_DIR}/parcels)
set(SUBMODELS ${CMAKE_CURRENT_SOURCE_DIR}/submodels)

# =====================================================================================
#                       Headers
# =====================================================================================
set(HEADERS_INTERMEDIATE
    # clouds
    ${CLOUDS}/baseClasses/carbonaceousCloud/carbonaceousCloud.H

    ${CLOUDS}/derived/basicCarbonaceousCloud/basicCarbonaceousCloud.H
    ${CLOUDS}/derived/basicReactingCloudSTFS/basicReactingCloudSTFS.H
    ${CLOUDS}/derived/basicReactingMultiphaseCloudSTFS/basicReactingMultiphaseCloudSTFS.H
    ${CLOUDS}/derived/basicThermoCloudSTFS/basicThermoCloudSTFS.H

    ${CLOUDS}/Templates/CarbonaceousCloud/CarbonaceousCloud.H
    ${CLOUDS}/Templates/ReactingCloudSTFS/ReactingCloudSTFS.H
    ${CLOUDS}/Templates/ReactingMultiphaseCloudSTFS/ReactingMultiphaseCloudSTFS.H
    ${CLOUDS}/Templates/ThermoCloudSTFS/ThermoCloudSTFS.H

    # parcels
    ${PARCELS}/derived/basicCarbonaceousParcel/basicCarbonaceousParcel.H
    ${PARCELS}/derived/basicReactingMultiphaseParcelSTFS/basicReactingMultiphaseParcelSTFS.H
    ${PARCELS}/derived/basicReactingParcelSTFS/basicReactingParcelSTFS.H
    ${PARCELS}/derived/basicThermoParcelSTFS/basicThermoParcelSTFS.H

    ${PARCELS}/include/makeCarbonaceousParcelDevolatilisationModels.H
    ${PARCELS}/include/makeCarbonaceousParcelHeatCapacityModels.H
    ${PARCELS}/include/makeCarbonaceousParcelSurfaceReactionModels.H
    ${PARCELS}/include/makeParcelSTFSHeatTransferModels.H

    ${PARCELS}/Templates/CarbonaceousParcel/CarbonaceousParcel.H
    ${PARCELS}/Templates/ReactingMultiphaseParcelSTFS/ReactingMultiphaseParcelSTFS.H
    ${PARCELS}/Templates/ReactingParcelSTFS/ReactingParcelSTFS.H
    ${PARCELS}/Templates/ReactingParcelSTFS/ReactingParcelSTFSTrackingDataI.H
    ${PARCELS}/Templates/ThermoParcelSTFS/ThermoParcelSTFS.H
    ${PARCELS}/Templates/ThermoParcelSTFS/ThermoParcelSTFSTrackingDataI.H

    # carbonaceous submodels
    ${SUBMODELS}/Carbonaceous/DevolatilizationModelSTFS/C2SM/C2SM.H
    ${SUBMODELS}/Carbonaceous/DevolatilizationModelSTFS/DevolatilisationModelSTFS/DevolatilisationModelSTFS.H
    ${SUBMODELS}/Carbonaceous/DevolatilizationModelSTFS/NoDevolatilisationSTFS/NoDevolatilisationSTFS.H
    ${SUBMODELS}/Carbonaceous/HeatCapacityModelSTFS/ConstantHeatCapacity/ConstantHeatCapacity.H
    ${SUBMODELS}/Carbonaceous/HeatCapacityModelSTFS/HeatCapacityFit/HeatCapacityFit.H
    ${SUBMODELS}/Carbonaceous/HeatCapacityModelSTFS/HeatCapacityModelSTFS/HeatCapacityModelSTFS.H
    ${SUBMODELS}/Carbonaceous/HeatCapacityModelSTFS/NoHeatCapacityModel/NoHeatCapacity.H
    ${SUBMODELS}/Carbonaceous/SurfaceReactionModelSTFS/COxidationKineticDiffusionLimitedRateSTFS/COxidationKineticDiffusionLimitedRateSTFS.H

    # cloud function objects
    ${SUBMODELS}/CloudFunctionObjects/SuperheatedIndex/SuperheatedIndex.H
    ${SUBMODELS}/CloudFunctionObjects/VoidMassFraction/VoidMassFraction.H
    ${SUBMODELS}/CloudFunctionObjects/ParticleVelocityField/ParticleVelocityField.H

    # cloud radiation submodels
    ${SUBMODELS}/CloudRadiation/wsggGreyCloudAbsorptionEmission/wsggGreyCloudAbsorptionEmission.H
)

# =====================================================================================
#                       Sources
# =====================================================================================
set(SOURCES_INTERMEDIATE
    # clouds
    ${CLOUDS}/baseClasses/carbonaceousCloud/carbonaceousCloud.C

    # derived parcels
    ${PARCELS}/derived/basicCarbonaceousParcel/defineBasicCarbonaceousParcel.C
    ${PARCELS}/derived/basicCarbonaceousParcel/makeBasicCarbonaceousParcelSubmodels.C
    ${PARCELS}/derived/basicReactingMultiphaseParcelSTFS/defineBasicReactingMultiphaseParcelSTFS.C
    ${PARCELS}/derived/basicReactingMultiphaseParcelSTFS/makeBasicReactingMultiphaseParcelSTFSSubmodels.C
    ${PARCELS}/derived/basicReactingParcelSTFS/defineBasicReactingParcelSTFS.C
    ${PARCELS}/derived/basicReactingParcelSTFS/makeBasicReactingParcelSTFSSubmodels.C
    ${PARCELS}/derived/basicThermoParcelSTFS/defineBasicThermoParcelSTFS.C
    ${PARCELS}/derived/basicThermoParcelSTFS/makeBasicThermoParcelSTFSSubmodels.C

    # cloud radiation submodels
    ${SUBMODELS}/CloudRadiation/wsggGreyCloudAbsorptionEmission/wsggGreyCloudAbsorptionEmission.C
)

add_library(lagrangianIntermediateSTFS SHARED
    ${HEADERS_FOAM}
    ${HEADERS_FOAM_LAGRANGIAN}
    ${HEADERS_INTERMEDIATE}
    ${SOURCES_INTERMEDIATE}
)

target_compile_options(lagrangianIntermediateSTFS PUBLIC -DNoRepository)

# =====================================================================================
#                       Link libraries
# =====================================================================================
target_link_libraries(lagrangianIntermediateSTFS PUBLIC
    OpenFOAM
    Pstream
    fileFormats
    finiteVolume
    fvOptions
    meshTools
    compressibleTurbulenceModels
    compressibleTransportModels
    distributionModels
    dynamicMesh
    dynamicFvMesh
    fluidThermophysicalModels
    incompressibleTransportModels
    lagrangian
    lagrangianIntermediate
    lagrangianTurbulence
    sampling
    SLGThermo
    specie
    surfMesh
    thermophysicalProperties
    topoChangerFvMesh
    turbulenceModels
    chemistryModel
    radiationModels
    ODE
    regionModels
    regionFaModels
    reactionThermophysicalModels
    surfaceFilmModels
    combustionModels
    finiteArea
    cantera_shared
    basicThermoSTFS
    canterathermo
    radiationModelsSTFS
    ${CMAKE_DL_LIBS}
)

# =====================================================================================
#                       Include directories
# =====================================================================================
target_include_directories(lagrangianIntermediateSTFS PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}

    # clouds
    ${CLOUDS}
    ${CLOUDS}/baseClasses
    ${CLOUDS}/baseClasses/carbonaceousCloud
    ${CLOUDS}/derived/basicCarbonaceousCloud
    ${CLOUDS}/derived/basicReactingCloudSTFS
    ${CLOUDS}/derived/basicReactingMultiphaseCloudSTFS
    ${CLOUDS}/derived/basicThermoCloudSTFS
    ${CLOUDS}/Templates
    ${CLOUDS}/Templates/CarbonaceousCloud
    ${CLOUDS}/Templates/ReactingCloudSTFS
    ${CLOUDS}/Templates/ReactingMultiphaseCloudSTFS
    ${CLOUDS}/Templates/ThermoCloudSTFS

    # parcels
    ${PARCELS}
    ${PARCELS}/derived
    ${PARCELS}/derived/basicCarbonaceousParcel
    ${PARCELS}/derived/basicReactingMultiphaseParcelSTFS
    ${PARCELS}/derived/basicReactingParcelSTFS
    ${PARCELS}/derived/basicThermoParcelSTFS
    ${PARCELS}/include
    ${PARCELS}/Templates
    ${PARCELS}/Templates/CarbonaceousParcel
    ${PARCELS}/Templates/ReactingMultiphaseParcelSTFS
    ${PARCELS}/Templates/ReactingParcelSTFS
    ${PARCELS}/Templates/ThermoParcelSTFS

    # submodels
    ${SUBMODELS}
    ${SUBMODELS}/Carbonaceous
    ${SUBMODELS}/Carbonaceous/DevolatilizationModelSTFS
    ${SUBMODELS}/Carbonaceous/DevolatilizationModelSTFS/C2SM
    ${SUBMODELS}/Carbonaceous/DevolatilizationModelSTFS/DevolatilisationModelSTFS
    ${SUBMODELS}/Carbonaceous/DevolatilizationModelSTFS/NoDevolatilisationSTFS
    ${SUBMODELS}/Carbonaceous/HeatCapacityModelSTFS
    ${SUBMODELS}/Carbonaceous/HeatCapacityModelSTFS/ConstantHeatCapacity
    ${SUBMODELS}/Carbonaceous/HeatCapacityModelSTFS/HeatCapacityFit
    ${SUBMODELS}/Carbonaceous/HeatCapacityModelSTFS/HeatCapacityModelSTFS
    ${SUBMODELS}/Carbonaceous/HeatCapacityModelSTFS/NoHeatCapacityModel
    ${SUBMODELS}/Carbonaceous/SurfaceReactionModelSTFS
    ${SUBMODELS}/Carbonaceous/SurfaceReactionModelSTFS/COxidationKineticDiffusionLimitedRateSTFS

    ${SUBMODELS}/CloudFunctionObjects
    ${SUBMODELS}/CloudFunctionObjects/SuperheatedIndex
    ${SUBMODELS}/CloudFunctionObjects/VoidMassFraction
    ${SUBMODELS}/CloudFunctionObjects/ParticleVelocityField

    ${SUBMODELS}/CloudRadiation
    ${SUBMODELS}/CloudRadiation/wsggGreyCloudAbsorptionEmission

    # CloudHandler.H only (header-only dependency: CloudHandler itself links
    # lagrangianIntermediateSTFS, so this cannot be a target_link_libraries
    # dependency without a cycle. The actual symbols resolve at the final
    # executable, which links both libraries.)
    ${CMAKE_SOURCE_DIR}/src/CloudHandler
)

install(TARGETS lagrangianIntermediateSTFS
    EXPORT lagrangianIntermediateSTFSConfig
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
install(FILES ${HEADERS_INTERMEDIATE} DESTINATION include)
