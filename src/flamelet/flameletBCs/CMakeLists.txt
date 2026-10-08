set(LIBNAME flameletBCs)

set(HEADERS_LOCAL
    ${CMAKE_CURRENT_SOURCE_DIR}/boundaryLookUp/boundaryLookUpFvPatchScalarField.H
)
set(SOURCES_LOCAL
    boundaryLookUp/boundaryLookUpFvPatchScalarField.C
)

add_library(${LIBNAME} SHARED
    ${HEADERS_FOAM}
    ${HEADERS_LOCAL}
    ${SOURCES_LOCAL}
)

target_link_libraries(${LIBNAME} PUBLIC
    finiteVolume
    ${CMAKE_DL_LIBS}
)

target_include_directories(${LIBNAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/boundaryLookUp
)

install(TARGETS ${LIBNAME}
    EXPORT ${LIBNAME}Config
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
install(FILES ${HEADERS_LOCAL} DESTINATION include)
