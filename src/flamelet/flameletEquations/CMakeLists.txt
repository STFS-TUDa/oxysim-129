set(LIBNAME flameletEquations)

set(HEADERS_LOCAL
    ${CMAKE_CURRENT_SOURCE_DIR}/flameletTypeDefs.H
    ${CMAKE_CURRENT_SOURCE_DIR}/inputVariableEquations.H
    ${CMAKE_CURRENT_SOURCE_DIR}/EquationHandler.H
)
set(SOURCES_LOCAL
    inputVariableEquations.C
    EquationHandler.C
)

add_library(${LIBNAME} SHARED
    ${HEADERS_FOAM}
    ${HEADERS_LOCAL}
    ${SOURCES_LOCAL}
)

target_compile_options(${LIBNAME} PUBLIC -DNoRepository)

target_link_libraries(${LIBNAME} PUBLIC
    turbulenceModels
    specie
    flameletThermo
    cloudHandler
    radiationModels
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
