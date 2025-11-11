#!/bin/bash
# benchmark_imports.sh
# Benchmark package import times for Apptainer vs uv/rix
#
# Measures:
# - NumPy import time (Python)
# - PyTorch import time (Python)
# - tidyverse load time (R)
# - Bioconductor package load times (R)
#
# Tests on both login nodes and compute nodes (via SLURM)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
RESULTS_FILE="${SCRIPT_DIR}/import_times.csv"

# Number of iterations for statistical significance
ITERATIONS=10

# Node type (login or compute)
NODE_TYPE="${1:-login}"

echo "=== Benchmark: Package Import Times ==="
echo "Node type: $NODE_TYPE"
echo "Iterations: $ITERATIONS"
echo "Results will be saved to: $RESULTS_FILE"
echo ""

# Initialize results file
if [ ! -f "$RESULTS_FILE" ]; then
    echo "method,node_type,package,iteration,time_ms" > "$RESULTS_FILE"
fi

# ============================================================================
# Helper Functions
# ============================================================================

# Measure import time in milliseconds
measure_python_import() {
    local method="$1"
    local package="$2"
    local cmd="$3"
    
    local start=$(date +%s%3N)
    eval "$cmd" > /dev/null 2>&1
    local end=$(date +%s%3N)
    echo $((end - start))
}

measure_r_load() {
    local method="$1"
    local package="$2"
    local cmd="$3"
    
    local start=$(date +%s%3N)
    eval "$cmd" > /dev/null 2>&1
    local end=$(date +%s%3N)
    echo $((end - start))
}

# Calculate statistics
calculate_stats() {
    local method="$1"
    local package="$2"
    local times=("${@:3}")
    
    # Calculate mean
    local sum=0
    for time in "${times[@]}"; do
        sum=$((sum + time))
    done
    local mean=$((sum / ${#times[@]}))
    
    # Calculate standard deviation
    local sq_diff_sum=0
    for time in "${times[@]}"; do
        local diff=$((time - mean))
        sq_diff_sum=$((sq_diff_sum + diff * diff))
    done
    local variance=$((sq_diff_sum / ${#times[@]}))
    local stddev=$(echo "scale=2; sqrt($variance)" | bc)
    
    # Find min and max
    local min=${times[0]}
    local max=${times[0]}
    for time in "${times[@]}"; do
        if [ $time -lt $min ]; then min=$time; fi
        if [ $time -gt $max ]; then max=$time; fi
    done
    
    printf "%-20s %-15s Mean: %6d ms  StdDev: %6.2f ms  Min: %6d ms  Max: %6d ms\n" \
        "$method" "$package" "$mean" "$stddev" "$min" "$max"
}

# ============================================================================
# Python Benchmarks
# ============================================================================

# Source init_env.sh
source "${ANALYSIS_DIR}/environments/init_env.sh" 2>/dev/null

echo "=== Python Package Imports ==="
echo ""

# ----------------------------------------------------------------------------
# NumPy Import - Apptainer
# ----------------------------------------------------------------------------
echo "Benchmarking NumPy import (Apptainer)..."
NUMPY_APPTAINER_TIMES=()

for i in $(seq 1 $ITERATIONS); do
    echo -n "  Iteration $i/$ITERATIONS... "
    
    time_taken=$(measure_python_import "apptainer" "numpy" \
        "sg grp_batch_effects -c 'apptainer exec --contain ~/groups/grp_batch_effects/remove-batch-effects.sif python -c \"import numpy\"'")
    
    NUMPY_APPTAINER_TIMES+=($time_taken)
    echo "$time_taken ms"
    echo "apptainer,$NODE_TYPE,numpy,$i,$time_taken" >> "$RESULTS_FILE"
done

echo ""
calculate_stats "Apptainer" "numpy" "${NUMPY_APPTAINER_TIMES[@]}"
echo ""

# ----------------------------------------------------------------------------
# NumPy Import - uv
# ----------------------------------------------------------------------------
echo "Benchmarking NumPy import (uv)..."
NUMPY_UV_TIMES=()

for i in $(seq 1 $ITERATIONS); do
    echo -n "  Iteration $i/$ITERATIONS... "
    
    time_taken=$(measure_python_import "uv" "numpy" \
        "bash -c 'source ${PYTHON_ENV}/bin/activate && python -c \"import numpy\"'")
    
    NUMPY_UV_TIMES+=($time_taken)
    echo "$time_taken ms"
    echo "uv,$NODE_TYPE,numpy,$i,$time_taken" >> "$RESULTS_FILE"
done

echo ""
calculate_stats "uv" "numpy" "${NUMPY_UV_TIMES[@]}"
echo ""

# ----------------------------------------------------------------------------
# PyTorch Import - Apptainer
# ----------------------------------------------------------------------------
echo "Benchmarking PyTorch import (Apptainer)..."
TORCH_APPTAINER_TIMES=()

for i in $(seq 1 $ITERATIONS); do
    echo -n "  Iteration $i/$ITERATIONS... "
    
    time_taken=$(measure_python_import "apptainer" "torch" \
        "sg grp_batch_effects -c 'apptainer exec --contain ~/groups/grp_batch_effects/remove-batch-effects.sif python -c \"import torch\"'")
    
    TORCH_APPTAINER_TIMES+=($time_taken)
    echo "$time_taken ms"
    echo "apptainer,$NODE_TYPE,torch,$i,$time_taken" >> "$RESULTS_FILE"
done

echo ""
calculate_stats "Apptainer" "torch" "${TORCH_APPTAINER_TIMES[@]}"
echo ""

# ----------------------------------------------------------------------------
# PyTorch Import - uv
# ----------------------------------------------------------------------------
echo "Benchmarking PyTorch import (uv)..."
TORCH_UV_TIMES=()

for i in $(seq 1 $ITERATIONS); do
    echo -n "  Iteration $i/$ITERATIONS... "
    
    time_taken=$(measure_python_import "uv" "torch" \
        "bash -c 'source ${PYTHON_ENV}/bin/activate && python -c \"import torch\"'")
    
    TORCH_UV_TIMES+=($time_taken)
    echo "$time_taken ms"
    echo "uv,$NODE_TYPE,torch,$i,$time_taken" >> "$RESULTS_FILE"
done

echo ""
calculate_stats "uv" "torch" "${TORCH_UV_TIMES[@]}"
echo ""

# ============================================================================
# R Benchmarks
# ============================================================================

NIX_ROOT="/grphome/grp_batch_effects/nix"
R_ENV_DIR="${ANALYSIS_DIR}/environments/r/batch-effects"
NIX_CACHE_OPTS="--option substituters 'https://cache.nixos.org https://rstats-on-nix.cachix.org' --option trusted-public-keys 'cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= rstats-on-nix.cachix.org-1:9cJb6nqYZgKqgH5XJQN8FPkXqKlGqKqJqKqKqKqKqKo='"

echo "=== R Package Loads ==="
echo ""

# Check if R environment exists
if [ ! -d "$R_ENV_DIR" ] || [ ! -f "$R_ENV_DIR/default.nix" ]; then
    echo "WARNING: R environment not found at $R_ENV_DIR"
    echo "Skipping R benchmarks. Run Phase 1 (authoring) first."
    echo ""
else
    # ----------------------------------------------------------------------------
    # tidyverse Load - Apptainer
    # ----------------------------------------------------------------------------
    echo "Benchmarking tidyverse load (Apptainer)..."
    TIDYVERSE_APPTAINER_TIMES=()
    
    for i in $(seq 1 $ITERATIONS); do
        echo -n "  Iteration $i/$ITERATIONS... "
        
        time_taken=$(measure_r_load "apptainer" "tidyverse" \
            "sg grp_batch_effects -c 'apptainer exec --contain ~/groups/grp_batch_effects/remove-batch-effects.sif Rscript -e \"library(tidyverse)\"'")
        
        TIDYVERSE_APPTAINER_TIMES+=($time_taken)
        echo "$time_taken ms"
        echo "apptainer,$NODE_TYPE,tidyverse,$i,$time_taken" >> "$RESULTS_FILE"
    done
    
    echo ""
    calculate_stats "Apptainer" "tidyverse" "${TIDYVERSE_APPTAINER_TIMES[@]}"
    echo ""
    
    # ----------------------------------------------------------------------------
    # tidyverse Load - rix
    # ----------------------------------------------------------------------------
    echo "Benchmarking tidyverse load (rix)..."
    TIDYVERSE_RIX_TIMES=()
    
    for i in $(seq 1 $ITERATIONS); do
        echo -n "  Iteration $i/$ITERATIONS... "
        
        time_taken=$(measure_r_load "rix" "tidyverse" \
            "$NIX_ROOT/nix-user-chroot $NIX_ROOT bash -c 'source ~/.nix-profile/etc/profile.d/nix.sh && cd $R_ENV_DIR && nix-shell $NIX_CACHE_OPTS --run \"Rscript -e \\\"library(tidyverse)\\\"\"'")
        
        TIDYVERSE_RIX_TIMES+=($time_taken)
        echo "$time_taken ms"
        echo "rix,$NODE_TYPE,tidyverse,$i,$time_taken" >> "$RESULTS_FILE"
    done
    
    echo ""
    calculate_stats "rix" "tidyverse" "${TIDYVERSE_RIX_TIMES[@]}"
    echo ""
    
    # ----------------------------------------------------------------------------
    # limma Load - Apptainer
    # ----------------------------------------------------------------------------
    echo "Benchmarking limma load (Apptainer)..."
    LIMMA_APPTAINER_TIMES=()
    
    for i in $(seq 1 $ITERATIONS); do
        echo -n "  Iteration $i/$ITERATIONS... "
        
        time_taken=$(measure_r_load "apptainer" "limma" \
            "sg grp_batch_effects -c 'apptainer exec --contain ~/groups/grp_batch_effects/remove-batch-effects.sif Rscript -e \"library(limma)\"'")
        
        LIMMA_APPTAINER_TIMES+=($time_taken)
        echo "$time_taken ms"
        echo "apptainer,$NODE_TYPE,limma,$i,$time_taken" >> "$RESULTS_FILE"
    done
    
    echo ""
    calculate_stats "Apptainer" "limma" "${LIMMA_APPTAINER_TIMES[@]}"
    echo ""
    
    # ----------------------------------------------------------------------------
    # limma Load - rix
    # ----------------------------------------------------------------------------
    echo "Benchmarking limma load (rix)..."
    LIMMA_RIX_TIMES=()
    
    for i in $(seq 1 $ITERATIONS); do
        echo -n "  Iteration $i/$ITERATIONS... "
        
        time_taken=$(measure_r_load "rix" "limma" \
            "$NIX_ROOT/nix-user-chroot $NIX_ROOT bash -c 'source ~/.nix-profile/etc/profile.d/nix.sh && cd $R_ENV_DIR && nix-shell $NIX_CACHE_OPTS --run \"Rscript -e \\\"library(limma)\\\"\"'")
        
        LIMMA_RIX_TIMES+=($time_taken)
        echo "$time_taken ms"
        echo "rix,$NODE_TYPE,limma,$i,$time_taken" >> "$RESULTS_FILE"
    done
    
    echo ""
    calculate_stats "rix" "limma" "${LIMMA_RIX_TIMES[@]}"
    echo ""
fi

# ============================================================================
# Summary
# ============================================================================
echo "=== Summary ==="
echo ""
echo "Results saved to: $RESULTS_FILE"
echo ""
echo "Comparison table:"
echo "--------------------------------------------------------------------------------"
printf "%-20s %-15s %10s %10s %10s %10s\n" "Method" "Package" "Mean (ms)" "StdDev" "Min (ms)" "Max (ms)"
echo "--------------------------------------------------------------------------------"

# Python packages
if [ ${#NUMPY_APPTAINER_TIMES[@]} -gt 0 ]; then
    calculate_stats "Apptainer" "numpy" "${NUMPY_APPTAINER_TIMES[@]}"
fi
if [ ${#NUMPY_UV_TIMES[@]} -gt 0 ]; then
    calculate_stats "uv" "numpy" "${NUMPY_UV_TIMES[@]}"
fi

if [ ${#TORCH_APPTAINER_TIMES[@]} -gt 0 ]; then
    calculate_stats "Apptainer" "torch" "${TORCH_APPTAINER_TIMES[@]}"
fi
if [ ${#TORCH_UV_TIMES[@]} -gt 0 ]; then
    calculate_stats "uv" "torch" "${TORCH_UV_TIMES[@]}"
fi

# R packages
if [ ${#TIDYVERSE_APPTAINER_TIMES[@]} -gt 0 ]; then
    calculate_stats "Apptainer" "tidyverse" "${TIDYVERSE_APPTAINER_TIMES[@]}"
fi
if [ ${#TIDYVERSE_RIX_TIMES[@]} -gt 0 ]; then
    calculate_stats "rix" "tidyverse" "${TIDYVERSE_RIX_TIMES[@]}"
fi

if [ ${#LIMMA_APPTAINER_TIMES[@]} -gt 0 ]; then
    calculate_stats "Apptainer" "limma" "${LIMMA_APPTAINER_TIMES[@]}"
fi
if [ ${#LIMMA_RIX_TIMES[@]} -gt 0 ]; then
    calculate_stats "rix" "limma" "${LIMMA_RIX_TIMES[@]}"
fi

echo "--------------------------------------------------------------------------------"
echo ""
echo "Done!"
