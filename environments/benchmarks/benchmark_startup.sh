#!/bin/bash
# benchmark_startup.sh
# Benchmark job startup times for Apptainer vs uv/rix
#
# Measures:
# - Apptainer startup (container + bind mounts)
# - uv activation time (Python-only)
# - nix-shell activation time (R-only)
# - Combined activation time (both environments)
#
# Runs each benchmark 10 times for statistical significance

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
RESULTS_FILE="${SCRIPT_DIR}/startup_times.csv"

# Number of iterations for statistical significance
ITERATIONS=10

echo "=== Benchmark: Job Startup Times ==="
echo "Iterations: $ITERATIONS"
echo "Results will be saved to: $RESULTS_FILE"
echo ""

# Initialize results file
echo "method,iteration,time_ms" > "$RESULTS_FILE"

# ============================================================================
# Helper Functions
# ============================================================================

# Measure time in milliseconds
time_ms() {
    local start=$(date +%s%3N)
    "$@" > /dev/null 2>&1
    local end=$(date +%s%3N)
    echo $((end - start))
}

# Calculate statistics
calculate_stats() {
    local method="$1"
    local times=("${@:2}")
    
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
    
    printf "%-30s Mean: %6d ms  StdDev: %6.2f ms  Min: %6d ms  Max: %6d ms\n" \
        "$method" "$mean" "$stddev" "$min" "$max"
}

# ============================================================================
# Benchmark 1: Apptainer Startup
# ============================================================================
echo "Benchmarking Apptainer startup..."
APPTAINER_TIMES=()

for i in $(seq 1 $ITERATIONS); do
    echo -n "  Iteration $i/$ITERATIONS... "
    
    # Measure time to start container and run simple command
    time_taken=$(time_ms sg grp_batch_effects -c \
        "apptainer exec --contain ~/groups/grp_batch_effects/remove-batch-effects.sif echo 'test'")
    
    APPTAINER_TIMES+=($time_taken)
    echo "$time_taken ms"
    echo "apptainer,$i,$time_taken" >> "$RESULTS_FILE"
done

echo ""
calculate_stats "Apptainer" "${APPTAINER_TIMES[@]}"
echo ""

# ============================================================================
# Benchmark 2: uv Activation (Python-only)
# ============================================================================
echo "Benchmarking uv activation (Python-only)..."
UV_TIMES=()

# Source init_env.sh once
source "${ANALYSIS_DIR}/environments/init_env.sh" 2>/dev/null

for i in $(seq 1 $ITERATIONS); do
    echo -n "  Iteration $i/$ITERATIONS... "
    
    # Measure time to activate Python environment and run simple command
    time_taken=$(time_ms bash -c "
        source '${PYTHON_ENV}/bin/activate' && \
        python -c 'print(\"test\")' \
    ")
    
    UV_TIMES+=($time_taken)
    echo "$time_taken ms"
    echo "uv_python,$i,$time_taken" >> "$RESULTS_FILE"
done

echo ""
calculate_stats "uv (Python-only)" "${UV_TIMES[@]}"
echo ""

# ============================================================================
# Benchmark 3: nix-shell Activation (R-only)
# ============================================================================
echo "Benchmarking nix-shell activation (R-only)..."
NIX_TIMES=()

NIX_ROOT="/grphome/grp_batch_effects/nix"
R_ENV_DIR="${ANALYSIS_DIR}/environments/r/batch-effects"
NIX_CACHE_OPTS="--option substituters 'https://cache.nixos.org https://rstats-on-nix.cachix.org' --option trusted-public-keys 'cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= rstats-on-nix.cachix.org-1:9cJb6nqYZgKqgH5XJQN8FPkXqKlGqKqJqKqKqKqKqKo='"

# Check if R environment exists
if [ ! -d "$R_ENV_DIR" ] || [ ! -f "$R_ENV_DIR/default.nix" ]; then
    echo "  WARNING: R environment not found at $R_ENV_DIR"
    echo "  Skipping R benchmarks. Run Phase 1 (authoring) first."
    echo ""
else
    for i in $(seq 1 $ITERATIONS); do
        echo -n "  Iteration $i/$ITERATIONS... "
        
        # Measure time to activate R environment and run simple command
        time_taken=$(time_ms bash -c "
            $NIX_ROOT/nix-user-chroot $NIX_ROOT bash -c \"
                source ~/.nix-profile/etc/profile.d/nix.sh && \
                cd '$R_ENV_DIR' && \
                nix-shell $NIX_CACHE_OPTS --run 'Rscript -e \\\"print(1)\\\"' \
            \" \
        ")
        
        NIX_TIMES+=($time_taken)
        echo "$time_taken ms"
        echo "nix_r,$i,$time_taken" >> "$RESULTS_FILE"
    done
    
    echo ""
    calculate_stats "nix-shell (R-only)" "${NIX_TIMES[@]}"
    echo ""
fi

# ============================================================================
# Benchmark 4: Combined Activation (Python + R)
# ============================================================================
echo "Benchmarking combined activation (Python + R)..."
COMBINED_TIMES=()

if [ ! -d "$R_ENV_DIR" ] || [ ! -f "$R_ENV_DIR/default.nix" ]; then
    echo "  WARNING: R environment not found. Skipping combined benchmarks."
    echo ""
else
    for i in $(seq 1 $ITERATIONS); do
        echo -n "  Iteration $i/$ITERATIONS... "
        
        # Measure time to activate both environments and run simple command
        time_taken=$(time_ms bash -c "
            source '${PYTHON_ENV}/bin/activate' && \
            $NIX_ROOT/nix-user-chroot $NIX_ROOT bash -c \"
                source ~/.nix-profile/etc/profile.d/nix.sh && \
                cd '$R_ENV_DIR' && \
                nix-shell $NIX_CACHE_OPTS --run 'python -c \\\"print(1)\\\" && Rscript -e \\\"print(1)\\\"' \
            \" \
        ")
        
        COMBINED_TIMES+=($time_taken)
        echo "$time_taken ms"
        echo "combined,$i,$time_taken" >> "$RESULTS_FILE"
    done
    
    echo ""
    calculate_stats "Combined (Python + R)" "${COMBINED_TIMES[@]}"
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
printf "%-30s %10s %10s %10s %10s\n" "Method" "Mean (ms)" "StdDev" "Min (ms)" "Max (ms)"
echo "--------------------------------------------------------------------------------"

if [ ${#APPTAINER_TIMES[@]} -gt 0 ]; then
    calculate_stats "Apptainer" "${APPTAINER_TIMES[@]}"
fi

if [ ${#UV_TIMES[@]} -gt 0 ]; then
    calculate_stats "uv (Python-only)" "${UV_TIMES[@]}"
fi

if [ ${#NIX_TIMES[@]} -gt 0 ]; then
    calculate_stats "nix-shell (R-only)" "${NIX_TIMES[@]}"
fi

if [ ${#COMBINED_TIMES[@]} -gt 0 ]; then
    calculate_stats "Combined (Python + R)" "${COMBINED_TIMES[@]}"
fi

echo "--------------------------------------------------------------------------------"
echo ""

# Calculate speedup if we have both Apptainer and uv times
if [ ${#APPTAINER_TIMES[@]} -gt 0 ] && [ ${#UV_TIMES[@]} -gt 0 ]; then
    apptainer_mean=$(echo "${APPTAINER_TIMES[@]}" | tr ' ' '\n' | awk '{sum+=$1} END {print sum/NR}')
    uv_mean=$(echo "${UV_TIMES[@]}" | tr ' ' '\n' | awk '{sum+=$1} END {print sum/NR}')
    speedup=$(echo "scale=2; $apptainer_mean / $uv_mean" | bc)
    echo "Speedup (Apptainer → uv): ${speedup}x faster"
    echo ""
fi

echo "Done!"
