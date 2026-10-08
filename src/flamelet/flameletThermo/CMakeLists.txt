set(LIBNAME flameletThermo)

set(HEADERS_LOCAL
    ${CMAKE_CURRENT_SOURCE_DIR}/flameletThermo.H
)
set(SOURCES_LOCAL
    flameletThermo.C
)

add_library(${LIBNAME} SHARED
    ${HEADERS_FOAM}
    ${HEADERS_LOCAL}
    ${SOURCES_LOCAL}
)

target_link_libraries(${LIBNAME} PUBLIC
    finiteVolume
    meshTools
    compressibleTransportModels
    fluidThermophysicalModels
    basicThermoSTFS
    flameletBCs
    flut
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
