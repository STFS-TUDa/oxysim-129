set(LIBNAME cloudHandler)

set(HEADERS_LOCAL
    ${CMAKE_CURRENT_SOURCE_DIR}/CloudHandler.H
)
set(SOURCES_LOCAL
    ${CMAKE_CURRENT_SOURCE_DIR}/CloudHandler.C
)

add_library(${LIBNAME} SHARED
    ${HEADERS_FOAM}
    ${HEADERS_FOAM_LAGRANGIAN}
    ${HEADERS_LOCAL}
    ${SOURCES_LOCAL}
)

target_compile_options(${LIBNAME} PUBLIC -DNoRepository)

target_link_libraries(${LIBNAME} PUBLIC
    lagrangianIntermediateSTFS
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
