#!/bin/bash
# conda_setup.sh
# Pure conda environment setup (no rv/uv)
# Usage: bash conda_setup.sh <project_name>

PROJECT_NAME="${1}"

if [ -z "$PROJECT_NAME" ]; then
    echo "Usage: bash conda_setup.sh <project_name>"
    echo ""
    echo "Available projects:"
    for dir in environments/*/; do
        if [ -f "$dir/environment.yml" ]; then
            basename "$dir"
        fi
    done
    exit 1
fi

PROJECT_DIR="environments/${PROJECT_NAME}"
ENV_FILE="${PROJECT_DIR}/environment.yml"

if [ ! -f "$ENV_FILE" ]; then
    echo "ERROR: environment.yml not found at $ENV_FILE"
    exit 1
fi

echo "=========================================="
echo "Pure Conda Environment Setup"
echo "=========================================="
echo "Project: $PROJECT_NAME"
echo "Environment file: $ENV_FILE"
echo ""

# Load conda if not available
if ! command -v conda &> /dev/null; then
    echo "Loading conda from modules..."
    for conda_module in miniconda3 anaconda3 conda; do
        if module -t avail "$conda_module" 2>&1 | grep -q "^${conda_module}"; then
            module load "$conda_module" 2>/dev/null || true
            if command -v conda &> /dev/null; then
                echo "✓ Loaded conda from module: $conda_module"
                break
            fi
        fi
    done
fi

if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not available"
    echo "Load conda manually: module load miniconda3"
    exit 1
fi

# Check if environment exists
if conda env list | grep -q "^${PROJECT_NAME} "; then
    echo "Conda environment '${PROJECT_NAME}' already exists"
    echo ""
    read -p "Recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        conda env remove -n "${PROJECT_NAME}" -y
    else
        echo "Keeping existing environment"
        echo ""
        echo "To activate:"
        echo "  conda activate ${PROJECT_NAME}"
        exit 0
    fi
fi

# Create environment
echo "Creating conda environment from $ENV_FILE..."
echo "This may take 10-30 minutes..."
echo ""

conda env create -f "$ENV_FILE"

echo ""
echo "=========================================="
echo "Installation Complete!"
echo "=========================================="
echo ""
echo "To activate:"
echo "  conda activate ${PROJECT_NAME}"
echo ""
echo "To deactivate:"
echo "  conda deactivate"
echo ""
echo "Environment includes:"
echo "  - R 4.4 + all R packages"
echo "  - Python 3.12 + scientific stack"
echo "  - System libraries (OpenBLAS, GLPK, etc.)"
echo ""
