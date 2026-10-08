set(LIBNAME radiationModelsSTFS)

set(RADIATIONMODEL ${CMAKE_CURRENT_SOURCE_DIR}/radiationModel)
set(SUBMODELS      ${CMAKE_CURRENT_SOURCE_DIR}/submodels)
set(PATCHFIELDS    ${CMAKE_CURRENT_SOURCE_DIR}/derivedFvPatchFields)

set(HEADERS_LOCAL
    ${RADIATIONMODEL}/radiationModelSTFS/radiationModelSTFS.H

    ${RADIATIONMODEL}/fvDOMwsgg/fvDOMwsgg/fvDOMwsgg.H
    ${RADIATIONMODEL}/fvDOMwsgg/fvDOMwsgg/fvDOMwsggI.H
    ${RADIATIONMODEL}/fvDOMwsgg/radiativeIntensityRaywsgg/radiativeIntensityRaywsgg.H
    ${RADIATIONMODEL}/fvDOMwsgg/radiativeIntensityRaywsgg/radiativeIntensityRaywsggI.H

    ${SUBMODELS}/wsggAbsorptionEmission/wsggAbsorptionEmission.H
    ${SUBMODELS}/wsggBordbarAbsorptionEmission/wsggBordbarAbsorptionEmission.H
    ${SUBMODELS}/wsggDorigonAbsorptionEmission/wsggDorigonAbsorptionEmission.H
    ${SUBMODELS}/wsggBinaryAbsorptionEmission/wsggBinaryAbsorptionEmission.H

    ${PATCHFIELDS}/wideBandDiffusiveRadiationwsgg/wideBandDiffusiveRadiationMixedFvPatchScalarFieldwsgg.H
)

set(SOURCES_LOCAL
    ${RADIATIONMODEL}/radiationModelSTFS/radiationModelSTFS.C

    ${RADIATIONMODEL}/fvDOMwsgg/fvDOMwsgg/fvDOMwsgg.C
    ${RADIATIONMODEL}/fvDOMwsgg/radiativeIntensityRaywsgg/radiativeIntensityRaywsgg.C

    ${SUBMODELS}/wsggAbsorptionEmission/wsggAbsorptionEmission.C
    ${SUBMODELS}/wsggBordbarAbsorptionEmission/wsggBordbarAbsorptionEmission.C
    ${SUBMODELS}/wsggDorigonAbsorptionEmission/wsggDorigonAbsorptionEmission.C
    ${SUBMODELS}/wsggBinaryAbsorptionEmission/wsggBinaryAbsorptionEmission.C

    ${PATCHFIELDS}/wideBandDiffusiveRadiationwsgg/wideBandDiffusiveRadiationMixedFvPatchScalarFieldwsgg.C
)

add_library(${LIBNAME} SHARED
    ${HEADERS_FOAM}
    ${HEADERS_LOCAL}
    ${SOURCES_LOCAL}
)

target_link_libraries(${LIBNAME} PUBLIC
    finiteVolume
    meshTools
    fluidThermophysicalModels
    radiationModels
    ${CMAKE_DL_LIBS}
)

target_include_directories(${LIBNAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}

    ${RADIATIONMODEL}/radiationModelSTFS
    ${RADIATIONMODEL}/fvDOMwsgg/fvDOMwsgg
    ${RADIATIONMODEL}/fvDOMwsgg/radiativeIntensityRaywsgg

    ${SUBMODELS}/wsggAbsorptionEmission
    ${SUBMODELS}/wsggBordbarAbsorptionEmission
    ${SUBMODELS}/wsggDorigonAbsorptionEmission
    ${SUBMODELS}/wsggBinaryAbsorptionEmission

    ${PATCHFIELDS}/wideBandDiffusiveRadiationwsgg
)

install(TARGETS ${LIBNAME}
    EXPORT ${LIBNAME}Config
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
install(FILES ${HEADERS_LOCAL} DESTINATION include)
