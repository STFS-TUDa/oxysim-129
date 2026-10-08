# =====================================================================================
#        Set name of the solver and directory! they have to be the same
# =====================================================================================
set(SOLVER directionalAverage)

# =====================================================================================
#        include directories to look for header files
# =====================================================================================
include_directories(
    ${UTLITIES}/${SOLVER} # important to add to solver dir
)

# =====================================================================================
#        include directories to look for libraries
# =====================================================================================
link_directories()

# =====================================================================================
#                       executable solver
# =====================================================================================
add_executable(${SOLVER}
    ${UTLITIES}/${SOLVER}/${SOLVER}.C
)

message("Added the following utility : ${SOLVER}.C")

# =====================================================================================
#        link required libraries
# =====================================================================================
target_link_libraries(${SOLVER} PUBLIC
    OpenFOAM
    pthread
    fileFormats
    finiteVolume
    fvOptions
    meshTools
    cloudHandler
    )

# =====================================================================================
#        add the solver to the install
# =====================================================================================
install(TARGETS ${SOLVER}
    ARCHIVE DESTINATION lib       # For static libraries
    LIBRARY DESTINATION lib       # For shared libraries
    RUNTIME DESTINATION bin       # For executables
)
