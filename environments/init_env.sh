#!/bin/bash
# Environment initialization script for uv/rix migration
# Sets up paths and environment variables for Python and R environments
#
# Usage:
#   source environments/init_env.sh              # Initialize with default settings
#   SHOW_STORAGE=1 source environments/init_env.sh  # Show storage usage
#   VERBOSE=1 source environments/init_env.sh    # Show detailed initialization

set -e

# ============================================================================
# Directory Structure (using explicit paths as per requirements)
# ============================================================================
# Note: /grphome is the actual path, ~/groups is a symlink
export SHARED_DIR="/grphome/grp_batch_effects"
export ANALYSIS_DIR="${HOME}/confounded_analysis"
export SCRIPTS_DIR="${ANALYSIS_DIR}/scripts"
export DATA_DIR="${SHARED_DIR}/data"
export OUTPUTS_DIR="${SHARED_DIR}/outputs"

# Scratch and archive storage (BYU RC storage tiers)
export SCRATCH_DIR="${SHARED_DIR}/nobackup/autodelete"  # 20 TiB, 12-week auto-delete
export ARCHIVE_DIR="${SHARED_DIR}/nobackup/archive"     # 20 TiB, long-term storage

# ============================================================================
# Python Environment Configuration (uv)
# ============================================================================
export PYTHON_ENV_SPEC="${ANALYSIS_DIR}/environments/python"
export PYTHON_ENV="${SHARED_DIR}/environments/python/.venv"  # Shared .venv in group storage
export UV_CACHE_DIR="${SHARED_DIR}/.uv_cache"                # Shared cache
export UV_PROJECT_ENVIRONMENT="${PYTHON_ENV}"                # Tell uv where to put .venv

# Add Python environment to PATH if it exists
if [ -d "${PYTHON_ENV}/bin" ]; then
    export PATH="${PYTHON_ENV}/bin:${PATH}"
fi

# Python environment variables for R argparse and reticulate
export PYTHON="${PYTHON_ENV}/bin/python"
export PYTHON3="${PYTHON_ENV}/bin/python"
export RETICULATE_PYTHON="${PYTHON_ENV}/bin/python"

# ============================================================================
# R Environment Configuration (rix/Nix)
# ============================================================================
export R_ENV_BATCH_EFFECTS="${ANALYSIS_DIR}/environments/r/batch-effects.nix"  # Main analysis (Bioc 3.21)
export R_ENV_COMBATSEQ="${ANALYSIS_DIR}/environments/r/combatseq.nix"          # ComBat-seq (Bioc 3.11)
export R_LIBS_USER="${SHARED_DIR}/environments/r/r-libs"                       # R packages in group storage

# Nix configuration (using rootless Nix from Task 1)
export NIX_WRAPPER="/grphome/grp_batch_effects/nix/nix-env"
export NIX_PATH="nixpkgs=channel:nixos-23.11"

# ============================================================================
# Directory Creation
# ============================================================================
# Create necessary directories if they don't exist (suppress errors if no permissions)

# Data directories
mkdir -p "${DATA_DIR}" "${DATA_DIR}/.cache" 2>/dev/null || true

# Output directories
mkdir -p "${OUTPUTS_DIR}/figures" "${OUTPUTS_DIR}/metrics" "${OUTPUTS_DIR}/tables" 2>/dev/null || true

# Environment directories in group storage
mkdir -p "${SHARED_DIR}/environments/python" 2>/dev/null || true
mkdir -p "${SHARED_DIR}/environments/r/r-libs" 2>/dev/null || true
mkdir -p "${UV_CACHE_DIR}" 2>/dev/null || true

# Scratch and archive directories (optional, for large temporary files)
mkdir -p "${SCRATCH_DIR}" 2>/dev/null || true
mkdir -p "${ARCHIVE_DIR}" 2>/dev/null || true

# ============================================================================
# Group Permissions Setup
# ============================================================================
# Set group-writable permissions for shared directories
# This allows all group members to use and update the environments
# Ignore errors if user doesn't own the directories

if [ "${VERBOSE:-0}" = "1" ]; then
    echo "Setting group permissions on shared directories..."
fi

# Set permissions on environment directories
chmod 775 "${SHARED_DIR}/environments" 2>/dev/null || true
chmod 775 "${SHARED_DIR}/environments/python" 2>/dev/null || true
chmod 775 "${SHARED_DIR}/environments/r" 2>/dev/null || true
chmod 775 "${SHARED_DIR}/environments/r/r-libs" 2>/dev/null || true
chmod 775 "${UV_CACHE_DIR}" 2>/dev/null || true

# Set setgid bit on directories so new files inherit group ownership
# This ensures all group members can access files created by others
chmod g+s "${SHARED_DIR}/environments" 2>/dev/null || true
chmod g+s "${SHARED_DIR}/environments/python" 2>/dev/null || true
chmod g+s "${SHARED_DIR}/environments/r" 2>/dev/null || true
chmod g+s "${UV_CACHE_DIR}" 2>/dev/null || true

# ============================================================================
# Storage Quota Checking Function
# ============================================================================
check_storage_quota() {
    echo ""
    echo "=== Storage Usage ==="
    echo ""
    
    # User Home (2 TiB quota, backed up)
    echo "User Home (2 TiB quota, backed up):"
    if [ -d "${ANALYSIS_DIR}" ]; then
        local home_size=$(du -sh "${ANALYSIS_DIR}" 2>/dev/null | cut -f1)
        echo "  ${ANALYSIS_DIR}: ${home_size}"
    else
        echo "  ${ANALYSIS_DIR}: Not found"
    fi
    
    # Group Home (2 TiB quota, backed up)
    echo ""
    echo "Group Home (2 TiB quota, backed up):"
    if [ -d "${SHARED_DIR}/environments" ]; then
        local env_size=$(du -sh "${SHARED_DIR}/environments" 2>/dev/null | cut -f1)
        echo "  Environments: ${env_size}"
    fi
    if [ -d "${UV_CACHE_DIR}" ]; then
        local cache_size=$(du -sh "${UV_CACHE_DIR}" 2>/dev/null | cut -f1)
        echo "  uv cache: ${cache_size}"
    fi
    if [ -d "/grphome/grp_batch_effects/nix" ]; then
        local nix_size=$(du -sh "/grphome/grp_batch_effects/nix" 2>/dev/null | cut -f1)
        echo "  Nix store: ${nix_size}"
    fi
    if [ -d "${SHARED_DIR}" ]; then
        local total_size=$(du -sh "${SHARED_DIR}" 2>/dev/null | cut -f1)
        echo "  Total group: ${total_size}"
    fi
    
    # Scratch (20 TiB quota, 12-week auto-delete)
    echo ""
    echo "Group Scratch (20 TiB quota, 12-week auto-delete):"
    if [ -d "${SCRATCH_DIR}" ]; then
        local scratch_size=$(du -sh "${SCRATCH_DIR}" 2>/dev/null | cut -f1)
        echo "  ${SCRATCH_DIR}: ${scratch_size}"
    else
        echo "  Not using scratch storage"
    fi
    
    # Archive (20 TiB quota, long-term storage)
    echo ""
    echo "Group Archive (20 TiB quota, long-term storage):"
    if [ -d "${ARCHIVE_DIR}" ]; then
        local archive_size=$(du -sh "${ARCHIVE_DIR}" 2>/dev/null | cut -f1)
        echo "  ${ARCHIVE_DIR}: ${archive_size}"
    else
        echo "  Not using archive storage"
    fi
    
    # Check for quota warnings (2 TiB = 2,000,000 MB)
    echo ""
    local home_usage_mb=$(du -sm "${ANALYSIS_DIR}" 2>/dev/null | cut -f1)
    if [ -n "$home_usage_mb" ] && [ "$home_usage_mb" -gt 1800000 ]; then
        echo "⚠️  WARNING: Home directory approaching 2 TiB quota (${home_usage_mb} MB used)"
    fi
    
    local group_usage_mb=$(du -sm "${SHARED_DIR}" 2>/dev/null | cut -f1)
    if [ -n "$group_usage_mb" ] && [ "$group_usage_mb" -gt 1800000 ]; then
        echo "⚠️  WARNING: Group directory approaching 2 TiB quota (${group_usage_mb} MB used)"
    fi
    
    echo ""
}

# ============================================================================
# Environment Verification Function
# ============================================================================
verify_environment() {
    local errors=0
    
    echo ""
    echo "=== Environment Verification ==="
    echo ""
    
    # Check Python environment
    if [ -d "${PYTHON_ENV}" ]; then
        echo "✓ Python environment exists: ${PYTHON_ENV}"
        if [ -f "${PYTHON_ENV}/bin/python" ]; then
            local py_version=$(${PYTHON_ENV}/bin/python --version 2>&1)
            echo "  Python version: ${py_version}"
        fi
    else
        echo "✗ Python environment not found: ${PYTHON_ENV}"
        echo "  Run: cd ${PYTHON_ENV_SPEC} && uv sync"
        errors=$((errors + 1))
    fi
    
    # Check R environment specs
    if [ -f "${R_ENV_BATCH_EFFECTS}" ]; then
        echo "✓ R batch-effects spec exists: ${R_ENV_BATCH_EFFECTS}"
    else
        echo "✗ R batch-effects spec not found: ${R_ENV_BATCH_EFFECTS}"
        errors=$((errors + 1))
    fi
    
    if [ -f "${R_ENV_COMBATSEQ}" ]; then
        echo "✓ R combatseq spec exists: ${R_ENV_COMBATSEQ}"
    else
        echo "⚠  R combatseq spec not found: ${R_ENV_COMBATSEQ} (optional)"
    fi
    
    # Check Nix wrapper
    if [ -f "${NIX_WRAPPER}" ]; then
        echo "✓ Nix wrapper exists: ${NIX_WRAPPER}"
    else
        echo "⚠  Nix wrapper not found: ${NIX_WRAPPER}"
        echo "  Nix environments may not be available"
    fi
    
    # Check directories
    if [ -d "${DATA_DIR}" ]; then
        echo "✓ Data directory exists: ${DATA_DIR}"
    else
        echo "✗ Data directory not found: ${DATA_DIR}"
        errors=$((errors + 1))
    fi
    
    if [ -d "${OUTPUTS_DIR}" ]; then
        echo "✓ Outputs directory exists: ${OUTPUTS_DIR}"
    else
        echo "✗ Outputs directory not found: ${OUTPUTS_DIR}"
        errors=$((errors + 1))
    fi
    
    echo ""
    if [ $errors -eq 0 ]; then
        echo "✓ All critical components verified"
    else
        echo "✗ Found ${errors} error(s) - some components missing"
    fi
    echo ""
    
    return $errors
}

# ============================================================================
# Display Initialization Status
# ============================================================================
if [ "${VERBOSE:-0}" = "1" ]; then
    echo ""
    echo "=== Environment Initialized ==="
    echo ""
    echo "Python Configuration:"
    echo "  Spec directory:     ${PYTHON_ENV_SPEC}"
    echo "  Virtual env:        ${PYTHON_ENV}"
    echo "  uv cache:           ${UV_CACHE_DIR}"
    echo ""
    echo "R Configuration:"
    echo "  Batch-effects env:  ${R_ENV_BATCH_EFFECTS}"
    echo "  ComBat-seq env:     ${R_ENV_COMBATSEQ}"
    echo "  R libraries:        ${R_LIBS_USER}"
    echo "  Nix wrapper:        ${NIX_WRAPPER}"
    echo ""
    echo "Data Directories:"
    echo "  Data:               ${DATA_DIR}"
    echo "  Outputs:            ${OUTPUTS_DIR}"
    echo "  Scratch:            ${SCRATCH_DIR}"
    echo "  Archive:            ${ARCHIVE_DIR}"
    echo ""
else
    echo "Environment initialized (use VERBOSE=1 for details)"
fi

# ============================================================================
# Optional Actions
# ============================================================================
# Show storage usage if requested
if [ "${SHOW_STORAGE:-0}" = "1" ]; then
    check_storage_quota
fi

# Verify environment if requested
if [ "${VERIFY:-0}" = "1" ]; then
    verify_environment
fi
