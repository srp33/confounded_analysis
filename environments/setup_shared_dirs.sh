#!/bin/bash
# Setup shared uv cache and .venv directories for group collaboration
# This script creates the necessary directories in /grphome with proper permissions

set -e

SHARED_DIR="/grphome/grp_batch_effects"
UV_CACHE_DIR="${SHARED_DIR}/.uv_cache"
PYTHON_ENV_DIR="${SHARED_DIR}/environments/python"

echo "Setting up shared Python environment directories..."
echo "Shared directory: ${SHARED_DIR}"

# Create directories
echo "Creating ${UV_CACHE_DIR}..."
mkdir -p "${UV_CACHE_DIR}"

echo "Creating ${PYTHON_ENV_DIR}..."
mkdir -p "${PYTHON_ENV_DIR}"

# Set group permissions (775 = rwxrwxr-x)
echo "Setting permissions to 775..."
chmod 775 "${UV_CACHE_DIR}" || echo "Warning: Could not set permissions on ${UV_CACHE_DIR}"
chmod 775 "${PYTHON_ENV_DIR}" || echo "Warning: Could not set permissions on ${PYTHON_ENV_DIR}"

# Set setgid bit so new files inherit group ownership
echo "Setting setgid bit for group inheritance..."
chmod g+s "${UV_CACHE_DIR}" || echo "Warning: Could not set setgid on ${UV_CACHE_DIR}"
chmod g+s "${PYTHON_ENV_DIR}" || echo "Warning: Could not set setgid on ${PYTHON_ENV_DIR}"

# Verify setup
echo ""
echo "Directory setup complete!"
echo ""
echo "Verification:"
ls -ld "${UV_CACHE_DIR}"
ls -ld "${PYTHON_ENV_DIR}"

echo ""
echo "Environment variables to set (add to init_env.sh):"
echo "export UV_CACHE_DIR=\"${UV_CACHE_DIR}\""
echo "export UV_PROJECT_ENVIRONMENT=\"${PYTHON_ENV_DIR}/.venv\""

echo ""
echo "Next steps:"
echo "1. Source init_env.sh to set environment variables"
echo "2. Run 'uv sync' from ${HOME}/confounded_analysis/environments/python/"
echo "3. The .venv will be created in ${PYTHON_ENV_DIR}/.venv/"
