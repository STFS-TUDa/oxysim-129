# CMake configuration for OpenFOAM-based builds.
# OpenFOAM's own bashrc/module must be sourced before running CMake —
# it sets WM_PROJECT_DIR, FOAM_SRC, FOAM_LIBBIN, WM_OPTIONS, FOAM_API, etc.
# Cantera can be provided via a CMake variable (preferred) or env var.

# ---------------------------------------------------------------------------
# Validate OpenFOAM environment
# ---------------------------------------------------------------------------
if(DEFINED ENV{WM_PROJECT_DIR})
    message(STATUS "OpenFOAM: " $ENV{WM_PROJECT_DIR})
else()
    message(FATAL_ERROR
        "OpenFOAM environment is not loaded.\n"
        "Source OpenFOAM's etc/bashrc (or load its module) before running CMake.")
endif()

set(OpenFOAM_VERSION $ENV{WM_PROJECT_VERSION})
set(OpenFOAM_DIR     $ENV{WM_PROJECT_DIR})
set(OpenFOAM_LIB_DIR $ENV{FOAM_LIBBIN})
set(OpenFOAM_SRC     $ENV{FOAM_SRC})

# ---------------------------------------------------------------------------
# Dependency root — CMake cache variable (preferred) with env var fallback.
# Pass on the command line: cmake -DCANTERA_ROOT=/opt/cantera
# ---------------------------------------------------------------------------
set(CANTERA_ROOT  "$ENV{CANTERA_ROOT}"  CACHE PATH "Cantera installation root")

if(CANTERA_ROOT)
    list(PREPEND CMAKE_PREFIX_PATH "${CANTERA_ROOT}")
    # Expose pkg-config files so pkg_check_modules(cantera) works without PKG_CONFIG_PATH
    set(ENV{PKG_CONFIG_PATH} "${CANTERA_ROOT}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    message(STATUS "CANTERA_ROOT: ${CANTERA_ROOT}")
else()
    message(WARNING "CANTERA_ROOT is not set. "
        "Pass -DCANTERA_ROOT=/path/to/cantera to cmake, or set the CANTERA_ROOT environment variable.")
endif()

# ---------------------------------------------------------------------------
# Build output paths
# Computed from WM_OPTIONS (set by OpenFOAM) + OXYSIM_INSTALL_DIR.
# Env vars FOAM_USER_APPBIN / FOAM_USER_LIBBIN override if set (backwards compat
# with existing source_*.sh scripts).
# ---------------------------------------------------------------------------
set(OXYSIM_INSTALL_DIR "default" CACHE STRING
    "Output sub-directory: build/plattforms-stfs-foam/<WM_OPTIONS>/<value>/")

set(FOAM_USER_APPBIN
    "${CMAKE_SOURCE_DIR}/build/plattforms-stfs-foam/$ENV{WM_OPTIONS}/${OXYSIM_INSTALL_DIR}/bin")
set(FOAM_USER_LIBBIN
    "${CMAKE_SOURCE_DIR}/build/plattforms-stfs-foam/$ENV{WM_OPTIONS}/${OXYSIM_INSTALL_DIR}/lib")

# ---------------------------------------------------------------------------
# Compiler flags (mirrors wmake -show-cxxflags for OpenFOAM ESI)
# ---------------------------------------------------------------------------
set(PATH_LIB_OPENMPI "openmpi-system")

add_compile_options(
    -std=c++20
    -DOPENFOAM=$ENV{FOAM_API}
    -pthread
    -DWM_ARCH_OPTION=64
    -DWM_DP
    -DWM_LABEL_SIZE=32
    -w
    -O3
    -DNoRepository
    -ftemplate-depth-100
    -fPIC)

option(SUPPRESS_WARNINGS "Suppress compiler warnings" ON)
if(SUPPRESS_WARNINGS)
    add_compile_options(
        -Wall
        -Wextra
        -Wold-style-cast
        -Wnon-virtual-dtor
        -Wno-unused-parameter
        -Wno-invalid-offsetof
        -Wno-attributes
        -Wno-unknown-pragmas)
endif()

add_definitions("${DEFINITIONS_COMPILE}")
if(APPLE)
    add_definitions("-Ddarwin64")
else()
    add_definitions("-Dlinux64")
endif()

# ---------------------------------------------------------------------------
# Global include paths
# ---------------------------------------------------------------------------
include_directories(PUBLIC
    ${CMAKE_SOURCE_DIR}/include

    # OpenFOAM source directories (alphabetical)
    ${OpenFOAM_SRC}/combustionModels/lnInclude
    ${OpenFOAM_SRC}/dynamicFvMesh/lnInclude
    ${OpenFOAM_SRC}/dynamicMesh/lnInclude
    ${OpenFOAM_SRC}/faOptions/lnInclude
    ${OpenFOAM_SRC}/fileFormats/lnInclude
    ${OpenFOAM_SRC}/finiteArea/lnInclude
    ${OpenFOAM_SRC}/finiteVolume/lnInclude
    ${OpenFOAM_SRC}/lagrangian/basic/lnInclude
    ${OpenFOAM_SRC}/lagrangian/distributionModels/lnInclude
    ${OpenFOAM_SRC}/lagrangian/intermediate/lnInclude
    ${OpenFOAM_SRC}/lagrangian/turbulence/lnInclude
    ${OpenFOAM_SRC}/lagrangian/spray/lnInclude
    ${OpenFOAM_SRC}/meshTools/lnInclude
    ${OpenFOAM_SRC}/ODE/lnInclude
    ${OpenFOAM_SRC}/OpenFOAM/lnInclude
    ${OpenFOAM_SRC}/OSspecific/POSIX/lnInclude
    ${OpenFOAM_SRC}/regionFaModels/lnInclude
    ${OpenFOAM_SRC}/regionModels/regionModel/lnInclude
    ${OpenFOAM_SRC}/regionModels/surfaceFilmModels/lnInclude
    ${OpenFOAM_SRC}/sampling/lnInclude
    ${OpenFOAM_SRC}/surfMesh/lnInclude
    ${OpenFOAM_SRC}/thermophysicalModels/basic/lnInclude
    ${OpenFOAM_SRC}/thermophysicalModels/radiation/lnInclude
    ${OpenFOAM_SRC}/thermophysicalModels/reactionThermo/lnInclude
    ${OpenFOAM_SRC}/thermophysicalModels/SLGThermo/lnInclude
    ${OpenFOAM_SRC}/thermophysicalModels/specie/lnInclude
    ${OpenFOAM_SRC}/thermophysicalModels/thermophysicalFunctions/lnInclude
    ${OpenFOAM_SRC}/thermophysicalModels/thermophysicalProperties/lnInclude
    ${OpenFOAM_SRC}/transportModels
    ${OpenFOAM_SRC}/transportModels/compressible/lnInclude
    ${OpenFOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude
    ${OpenFOAM_SRC}/TurbulenceModels/incompressible/lnInclude
    ${OpenFOAM_SRC}/TurbulenceModels/compressible/lnInclude

    # Third-party includes (resolved via cache variable above)
    ${CANTERA_ROOT}/include
)

# ---------------------------------------------------------------------------
# Global link paths
# ---------------------------------------------------------------------------
link_directories(
    ${OpenFOAM_LIB_DIR}
    ${OpenFOAM_LIB_DIR}/dummy
    ${OpenFOAM_LIB_DIR}/${PATH_LIB_OPENMPI}
    ${FOAM_USER_LIBBIN}
    ${CANTERA_ROOT}/lib
)
