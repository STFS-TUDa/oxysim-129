set(LIBNAME basicThermoSTFS)

set(HEADERS_LOCAL
    ${CMAKE_CURRENT_SOURCE_DIR}/basicThermoSTFS.H
)
set(SOURCES_LOCAL
    ${CMAKE_CURRENT_SOURCE_DIR}/basicThermoSTFS.C
)

add_library(${LIBNAME} SHARED
    ${HEADERS_FOAM}
    ${HEADERS_LOCAL}
    ${SOURCES_LOCAL}
)

target_link_libraries(${LIBNAME} PUBLIC
    OpenFOAM
    Pstream
    compressibleTransportModels
    fluidThermophysicalModels
    specie
    finiteVolume
    ${CMAKE_DL_LIBS}
)

target_include_directories(${LIBNAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

install(TARGETS ${LIBNAME}
    EXPORT ${LIBNAME}Config
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
install(FILES ${HEADERS_LOCAL} DESTINATION include)
