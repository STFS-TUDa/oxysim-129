set(SOLVER flameletFoam)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(${SOLVER}
    ${CMAKE_CURRENT_SOURCE_DIR}/${SOLVER}.C
)

message("Added the following solver: ${SOLVER}.C")

target_link_libraries(${SOLVER} PUBLIC
    Pstream
    fvOptions
    fileFormats
    finiteVolume
    compressibleTurbulenceModels
    compressibleTransportModels
    distributionModels
    dynamicMesh
    dynamicFvMesh
    fluidThermophysicalModels
    incompressibleTransportModels
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
    basicThermoSTFS
    flameletThermo
    flameletEquations
    basicBenchmarkUtility
)

install(TARGETS ${SOLVER}
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
