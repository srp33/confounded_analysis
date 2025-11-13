#!/bin/bash
# create_conda_r_env.sh
# Create conda environments with R + system libraries for use with rv
# Usage: bash create_conda_r_env.sh <r_version>

set -e

R_VERSION="${1}"

if [ -z "$R_VERSION" ]; then
    echo "Usage: bash create_conda_r_env.sh <r_version>"
    echo ""
    echo "Examples:"
    echo "  bash create_conda_r_env.sh 4.0"
    echo "  bash create_conda_r_env.sh 4.4"
    echo "  bash create_conda_r_env.sh 4.5"
    exit 1
fi

# Normalize version (4.0 -> 4.0, 4.4.0 -> 4.4)
R_VERSION_SHORT=$(echo "$R_VERSION" | cut -d. -f1,2)
CONDA_ENV_NAME="rv-r${R_VERSION_SHORT}-syslibs"

echo "=========================================="
echo "Creating Conda R Environment"
echo "=========================================="
echo "R Version: ${R_VERSION}"
echo "Conda env name: ${CONDA_ENV_NAME}"
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found"
    echo ""
    echo "Load conda first:"
    echo "  module load anaconda3"
    echo "Or:"
    echo "  module load miniconda3"
    exit 1
fi

# Check if environment already exists
if conda env list | grep -q "^${CONDA_ENV_NAME} "; then
    echo "Conda environment '${CONDA_ENV_NAME}' already exists"
    echo ""
    read -p "Do you want to recreate it? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Keeping existing environment"
        exit 0
    fi
    echo "Removing existing environment..."
    conda env remove -n "${CONDA_ENV_NAME}" -y
fi

echo "Creating conda environment with R ${R_VERSION} and system libraries..."
echo ""

# Create environment with R and essential system libraries
# Note: FlexiBLAS not available in conda, using OpenBLAS instead
echo "Searching for R ${R_VERSION} packages..."
conda create -n "${CONDA_ENV_NAME}" -y \
    -c conda-forge \
    "r-base=${R_VERSION}.*" \
    openblas \
    glpk \
    libcurl \
    libxml2 \
    openssl \
    gfortran_linux-64 \
    make \
    cmake

echo ""
echo "=========================================="
echo "Installation Complete!"
echo "=========================================="
echo ""
echo "Conda environment '${CONDA_ENV_NAME}' created with:"
echo "  - R ${R_VERSION}"
echo "  - OpenBLAS (linear algebra - replaces FlexiBLAS)"
echo "  - GLPK (optimization)"
echo "  - libcurl, libxml2, openssl (networking)"
echo "  - Build tools (gfortran, make, cmake)"
echo ""
echo "This environment provides R and system libraries."
echo "R packages will be managed by rv (fast, up-to-date)."
echo ""
echo "Next steps:"
echo "  1. Update rproject.toml to use r_version = \"${R_VERSION_SHORT}\""
echo "  2. Activate: source environments/load_envs.sh <project_name>"
echo ""
