# =====================================================================================
#                       Subdirectory — build order follows dependency graph
# =====================================================================================

# Thermo foundation
add_subdirectory(thermo/basicThermoSTFS)
add_subdirectory(thermo/canteraThermo)

# Radiation (must precede lagrangian, which links it for cloud absorption/emission)
add_subdirectory(radiation)

# Lagrangian (must precede flameletEquations which links cloudHandler)
add_subdirectory(lagrangian)
add_subdirectory(CloudHandler)

# Flamelet group
add_subdirectory(flamelet/flameletBCs)
add_subdirectory(flamelet/flameletThermo)
add_subdirectory(flamelet/flameletEquations)

# Utilities
add_subdirectory(benchmarkUtility)

