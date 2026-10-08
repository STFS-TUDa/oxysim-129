# use "cmake --preset ninja ../." to configure this project
cmake_minimum_required(VERSION 3.28)
# Define the project name and version
project(OxySim-129 VERSION 1.0 LANGUAGES CXX C)

# Do not build in the source directory! bad practice, create a build directory with "mkdir build"
if (CMAKE_BINARY_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message(FATAL_ERROR "Building in-source is not supported! Create a build dir and remove ${CMAKE_SOURCE_DIR}/CMakeCache.txt")
endif()

# =====================================================================================
#             include cmake configuration files
# =====================================================================================
# Set the CMAKE_MODULE_PATH to include the custom 'cmake' directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# include cmake dependent option
include(CMakeDependentOption)

include(customMacros)
include(OpenFOAM)

set(CMAKE_CXX_STANDARD_REQUIRED True)

# Use Gold as the linker ->  faster than ld
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=gold")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=gold")

# set output directories based on OpenFOAM
message(STATUS "${FOAM_USER_APPBIN}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${FOAM_USER_APPBIN})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${FOAM_USER_LIBBIN})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${FOAM_USER_LIBBIN})

# install output
# install output
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/../install CACHE PATH "Default install path" FORCE)
endif()


get_env_var(COMPILE_LAGRANGIAN_WITH_TESTS)

# set RUNPATH so that the binaries are not dependent on LD_LIBRARY_PATH
# to be set up correctly. The binaries will always search for shared libraries
# in the directories they were linked against during build time.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# ensure catch2 is fetched before any other subproject is included
if(COMPILE_LAGRANGIAN_WITH_TESTS)
    message(STATUS "Compile with Catch2")
    include(FetchContent)
    FetchContent_Declare(
        catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.7.0
    )
    FetchContent_MakeAvailable(catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/src)
    include(CTest)
    include(Catch)
endif()

# =====================================================================================
#                       Add flut-reader (builds the flut library)
# =====================================================================================
add_subdirectory(extern/flut-reader)

# =====================================================================================
#                       Add libraries
# =====================================================================================
set(LIBRARIES ${PROJECT_SOURCE_DIR}/src)
add_subdirectory(${LIBRARIES})

# =====================================================================================
#                       Add applications
# =====================================================================================
set(APPLICATIONS ${PROJECT_SOURCE_DIR}/applications)
add_subdirectory(${APPLICATIONS})

# =====================================================================================
#                       Add tests
# =====================================================================================
if(COMPILE_LAGRANGIAN_WITH_TESTS)
    add_subdirectory(tests)
endif()

include_directories(${PROJECT_SOURCE_DIR}/src)

# =====================================================================================
#                       Install scripts
# =====================================================================================
install(DIRECTORY ${CMAKE_SOURCE_DIR}/scripts DESTINATION ${CMAKE_INSTALL_PREFIX})
