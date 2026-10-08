# =====================================================================================
#        Set name of the solver and directory! they have to be the same
# =====================================================================================
set(SOLVER flameletCloudFoam)

# =====================================================================================
#        include directories to look for header files 
# =====================================================================================
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(${SOLVER}
    ${CMAKE_CURRENT_SOURCE_DIR}/${SOLVER}.C
)

message("Added the following solver : ${SOLVER}.C")

# =====================================================================================
#        link required libraries
# =====================================================================================
target_link_libraries(${SOLVER} PUBLIC
    pthread
    cloudHandler
    flameletThermo
    flameletBCs
    flameletEquations
    flut
    radiationModels
    radiationModelsSTFS
    cantera_shared
)

# =====================================================================================
#        add the solver to the install
# =====================================================================================
install(TARGETS ${SOLVER}
    ARCHIVE DESTINATION lib       # For static libraries
    LIBRARY DESTINATION lib       # For shared libraries
    RUNTIME DESTINATION bin       # For executables 
)
