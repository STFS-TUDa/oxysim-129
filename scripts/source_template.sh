#!/bin/bash
# Template environment-setup script for OxySim-129.
#
# Copy this file (e.g. `cp scripts/source_template.sh scripts/source_local.sh`),
# fill in / adjust what applies to your system below, then source it
# before building or running OxySim-129:
#
#   source scripts/source_local.sh
#
# See the top-level README's "Installation"/"Running" sections for more
# detail on each of these paths.

# --- Fill these in for your system --------------------------------------

# OpenFOAM installation (its etc/bashrc)
OPENFOAM_BASHRC="/path/to/OpenFOAM-v2512/etc/bashrc"

# Cantera: loaded via a module system by default (edit the module name for
# your system). If it's not available as a module, comment this line out and
# uncomment/fill in the export line below instead, pointing at its install
# prefix (containing include/ and lib/).
module load cantera/<version>
# export CANTERA_ROOT="/path/to/cantera"

# Sub-directory name passed to -DOXYSIM_INSTALL_DIR at build time, if any
# (leave as "default" if you didn't set it)
OXYSIM_INSTALL_DIR="default"

# --------------------------------------------------------------------------

source "${OPENFOAM_BASHRC}"

# Repository root = one level up from wherever this script lives.
# Exported because tests/integration/*/foam_template/Allrun scripts source
# $OXYSIM_REPO_ROOT/bin/tools/STFSRunFunctions - required even though pytest
# copies foam_template/ into a temp/ directory that doesn't include bin/tools/
export OXYSIM_REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

export FOAM_USER_APPBIN="${OXYSIM_REPO_ROOT}/build/plattforms-stfs-foam/${WM_OPTIONS}/${OXYSIM_INSTALL_DIR}/bin"
export FOAM_USER_LIBBIN="${OXYSIM_REPO_ROOT}/build/plattforms-stfs-foam/${WM_OPTIONS}/${OXYSIM_INSTALL_DIR}/lib"
export PATH="${FOAM_USER_APPBIN}:${PATH}"
export LD_LIBRARY_PATH="${FOAM_USER_LIBBIN}:${CANTERA_ROOT}/lib:${LD_LIBRARY_PATH}"

echo "OxySim-129 environment loaded (repo: ${OXYSIM_REPO_ROOT}, WM_OPTIONS: ${WM_OPTIONS})"
