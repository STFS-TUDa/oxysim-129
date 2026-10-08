set(LIBNAME canterathermo)

set(HEADERS_LOCAL
    ${CMAKE_CURRENT_SOURCE_DIR}/canteraSpecieMixture.H
    ${CMAKE_CURRENT_SOURCE_DIR}/canteraThermo.H
)
set(SOURCES_LOCAL
    canteraSpecieMixture.C
    canteraThermo.C
)

add_library(${LIBNAME} SHARED
    ${HEADERS_FOAM}
    ${HEADERS_LOCAL}
    ${SOURCES_LOCAL}
)

find_package(PkgConfig REQUIRED)
pkg_check_modules(CANTERA REQUIRED cantera)

target_link_libraries(${LIBNAME} PUBLIC
    finiteVolume
    meshTools
    basicThermoSTFS
    ${CANTERA_LIBRARIES}
    ${CMAKE_DL_LIBS}
)

target_link_directories(${LIBNAME} PUBLIC
    ${CANTERA_LIBRARY_DIRS}
)

target_include_directories(${LIBNAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CANTERA_INCLUDEDIR}
)

install(TARGETS ${LIBNAME}
    EXPORT ${LIBNAME}Config
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
install(FILES ${HEADERS_LOCAL} DESTINATION include)
