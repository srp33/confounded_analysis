#!/bin/bash
# benchmark_nix_build_optimization.sh
# Benchmark nix-build optimization performance
#
# Measures:
# - Activation time with default.nix (baseline: ~300s)
# - Activation time with ./result (optimized: target <5s)
# - Package loading time with optimized activation
# - Speedup achieved
#
# Runs 10 iterations for statistical significance

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
RESULTS_FILE="${SCRIPT_DIR}/nix_build_optimization_results.csv"

# Number of iterations
ITERATIONS=10

# Nix configuration
NIX_ROOT="/grphome/grp_batch_effects/nix"
NIX_CHROOT_CMD="$NIX_ROOT/nix-user-chroot $NIX_ROOT"
NIX_CACHE_OPTS="--option substituters 'https://cache.nixos.org https://rstats-on-nix.cachix.org' --option trusted-public-keys 'cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= rstats-on-nix.cachix.org-1:9cJb6nqYZgKqgH5XJQN8FPkXqKlGqKqJqKqKqKqKqKo='"

# R environment directories
BATCH_EFFECTS_DIR="${ANALYSIS_DIR}/environments/r/batch-effects"
COMBATSEQ_DIR="${ANALYSIS_DIR}/environments/r/combatseq"

echo "============================================"
echo "Benchmark: nix-build Optimization"
echo "============================================"
echo "Iterations: $ITERATIONS"
echo "Results: $RESULTS_FILE"
echo ""

# Initialize results file
echo "environment,method,iteration,time_ms" > "$RESULTS_FILE"

# Initialize arrays
RESULT_TIMES=()
PACKAGE_TIMES=()
COMBATSEQ_TIMES=()
DEFAULT_TIMES=()

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
    local label="$1"
    shift
    local times=("$@")
    
    if [ ${#times[@]} -eq 0 ]; then
        echo "No data for $label"
        return
    fi
    
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
    local stddev=$(echo "scale=2; sqrt($variance)" | bc 2>/dev/null || echo "0")
    
    # Find min and max
    local min=${times[0]}
    local max=${times[0]}
    for time in "${times[@]}"; do
        if [ $time -lt $min ]; then min=$time; fi
        if [ $time -gt $max ]; then max=$time; fi
    done
    
    printf "%-40s Mean: %8d ms  StdDev: %8.2f ms  Min: %8d ms  Max: %8d ms\n" \
        "$label" "$mean" "$stddev" "$min" "$max"
}

# Build environment if result doesn't exist
ensure_result_exists() {
    local env_dir="$1"
    local env_name=$(basename "$env_dir")
    
    if [ ! -L "$env_dir/result" ]; then
        echo "⚠ WARNING: No ./result symlink found for $env_name"
        echo ""
        echo "The R environment needs to be built first using:"
        echo "  cd $env_dir"
        echo "  ../../build_nix_env.sh"
        echo ""
        echo "This is a one-time operation that takes ~5-10 minutes."
        echo "Once complete, re-run this benchmark script."
        echo ""
        return 1
    else
        echo "✓ Found existing result for $env_name"
        echo "  Target: $(readlink "$env_dir/result")"
        echo ""
    fi
}

# ============================================================================
# Benchmark batch-effects environment
# ============================================================================

if [ -d "$BATCH_EFFECTS_DIR" ] && [ -f "$BATCH_EFFECTS_DIR/default.nix" ]; then
    echo "============================================"
    echo "Benchmarking: batch-effects environment"
    echo "============================================"
    echo ""
    
    # Ensure result exists
    ensure_result_exists "$BATCH_EFFECTS_DIR" || {
        echo "Skipping batch-effects benchmarks"
        echo ""
    }
    
    if [ -L "$BATCH_EFFECTS_DIR/result" ]; then
        # Benchmark 1: Activation with result (optimized)
        echo "Test 1: Activation with ./result (optimized)"
        RESULT_TIMES=()
        
        for i in $(seq 1 $ITERATIONS); do
            echo -n "  Iteration $i/$ITERATIONS... "
            
            time_taken=$(time_ms bash -c "
                $NIX_CHROOT_CMD bash -c \"
                    source ~/.nix-profile/etc/profile.d/nix.sh && \
                    cd '$BATCH_EFFECTS_DIR' && \
                    nix-shell result $NIX_CACHE_OPTS --run 'R --version'
                \"
            ")
            
            RESULT_TIMES+=($time_taken)
            echo "$time_taken ms"
            echo "batch-effects,result,$i,$time_taken" >> "$RESULTS_FILE"
        done
        
        echo ""
        calculate_stats "batch-effects (result)" "${RESULT_TIMES[@]}"
        echo ""
        
        # Benchmark 2: Package loading with result
        echo "Test 2: Package loading (tidyverse)"
        PACKAGE_TIMES=()
        
        for i in $(seq 1 $ITERATIONS); do
            echo -n "  Iteration $i/$ITERATIONS... "
            
            time_taken=$(time_ms bash -c "
                $NIX_CHROOT_CMD bash -c \"
                    source ~/.nix-profile/etc/profile.d/nix.sh && \
                    cd '$BATCH_EFFECTS_DIR' && \
                    nix-shell result $NIX_CACHE_OPTS --run 'Rscript -e \\\"library(tidyverse)\\\"'
                \"
            ")
            
            PACKAGE_TIMES+=($time_taken)
            echo "$time_taken ms"
            echo "batch-effects,package_load,$i,$time_taken" >> "$RESULTS_FILE"
        done
        
        echo ""
        calculate_stats "batch-effects (package load)" "${PACKAGE_TIMES[@]}"
        echo ""
        
        # Benchmark 3: Baseline with default.nix (optional, slow)
        echo "Test 3: Baseline with default.nix (optional)"
        read -p "Run baseline test? Takes ~5 minutes per iteration. [y/N] " -n 1 -r
        echo
        
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            echo "Running baseline test (this will take ~50 minutes for 10 iterations)..."
            echo "Press Ctrl+C to skip after a few iterations if needed."
            echo ""
            
            DEFAULT_TIMES=()
            
            for i in $(seq 1 $ITERATIONS); do
                echo -n "  Iteration $i/$ITERATIONS... "
                
                time_taken=$(time_ms bash -c "
                    $NIX_CHROOT_CMD bash -c \"
                        source ~/.nix-profile/etc/profile.d/nix.sh && \
                        cd '$BATCH_EFFECTS_DIR' && \
                        nix-shell default.nix $NIX_CACHE_OPTS --run 'R --version'
                    \"
                ")
                
                DEFAULT_TIMES+=($time_taken)
                echo "$time_taken ms"
                echo "batch-effects,default,$i,$time_taken" >> "$RESULTS_FILE"
            done
            
            echo ""
            calculate_stats "batch-effects (default.nix)" "${DEFAULT_TIMES[@]}"
            echo ""
            
            # Calculate speedup
            if [ ${#DEFAULT_TIMES[@]} -gt 0 ] && [ ${#RESULT_TIMES[@]} -gt 0 ]; then
                default_mean=$(echo "${DEFAULT_TIMES[@]}" | tr ' ' '\n' | awk '{sum+=$1} END {print sum/NR}')
                result_mean=$(echo "${RESULT_TIMES[@]}" | tr ' ' '\n' | awk '{sum+=$1} END {print sum/NR}')
                speedup=$(echo "scale=2; $default_mean / $result_mean" | bc)
                
                echo "Speedup: ${speedup}x faster (default.nix → result)"
                echo ""
            fi
        else
            echo "Skipped baseline test."
            echo ""
        fi
    fi
else
    echo "WARNING: batch-effects environment not found at $BATCH_EFFECTS_DIR"
    echo "Skipping batch-effects benchmarks."
    echo ""
fi

# ============================================================================
# Benchmark combatseq environment
# ============================================================================

if [ -d "$COMBATSEQ_DIR" ] && [ -f "$COMBATSEQ_DIR/default.nix" ]; then
    echo "============================================"
    echo "Benchmarking: combatseq environment"
    echo "============================================"
    echo ""
    
    # Ensure result exists
    ensure_result_exists "$COMBATSEQ_DIR" || {
        echo "Skipping combatseq benchmarks"
        echo ""
    }
    
    if [ -L "$COMBATSEQ_DIR/result" ]; then
        # Benchmark: Activation with result
        echo "Test: Activation with ./result (optimized)"
        COMBATSEQ_TIMES=()
        
        for i in $(seq 1 $ITERATIONS); do
            echo -n "  Iteration $i/$ITERATIONS... "
            
            time_taken=$(time_ms bash -c "
                $NIX_CHROOT_CMD bash -c \"
                    source ~/.nix-profile/etc/profile.d/nix.sh && \
                    cd '$COMBATSEQ_DIR' && \
                    nix-shell result $NIX_CACHE_OPTS --run 'R --version'
                \"
            ")
            
            COMBATSEQ_TIMES+=($time_taken)
            echo "$time_taken ms"
            echo "combatseq,result,$i,$time_taken" >> "$RESULTS_FILE"
        done
        
        echo ""
        calculate_stats "combatseq (result)" "${COMBATSEQ_TIMES[@]}"
        echo ""
    fi
else
    echo "WARNING: combatseq environment not found at $COMBATSEQ_DIR"
    echo "Skipping combatseq benchmarks."
    echo ""
fi

# ============================================================================
# Summary
# ============================================================================

echo "============================================"
echo "Summary"
echo "============================================"
echo ""
echo "Results saved to: $RESULTS_FILE"
echo ""

# Check success criteria
if [ "${#RESULT_TIMES[@]}" -gt 0 ] 2>/dev/null; then
    result_mean=$(echo "${RESULT_TIMES[@]}" | tr ' ' '\n' | awk '{sum+=$1} END {print sum/NR}')
    result_mean_s=$(echo "scale=2; $result_mean / 1000" | bc)
    
    echo "batch-effects activation: ${result_mean_s}s (mean)"
    
    if [ $(echo "$result_mean < 5000" | bc) -eq 1 ]; then
        echo "✓ PASS: Activation time < 5s (target achieved)"
    else
        echo "⚠ WARNING: Activation time > 5s (expected ~2s)"
    fi
    echo ""
else
    echo "⚠ No batch-effects activation data collected"
    echo ""
fi

if [ "${#PACKAGE_TIMES[@]}" -gt 0 ] 2>/dev/null; then
    package_mean=$(echo "${PACKAGE_TIMES[@]}" | tr ' ' '\n' | awk '{sum+=$1} END {print sum/NR}')
    package_mean_s=$(echo "scale=2; $package_mean / 1000" | bc)
    
    echo "Package loading (tidyverse): ${package_mean_s}s (mean)"
    
    if [ $(echo "$package_mean < 10000" | bc) -eq 1 ]; then
        echo "✓ PASS: Package loading < 10s"
    else
        echo "⚠ WARNING: Package loading > 10s"
    fi
    echo ""
else
    echo "⚠ No package loading data collected"
    echo ""
fi

if [ "${#COMBATSEQ_TIMES[@]}" -gt 0 ] 2>/dev/null; then
    combatseq_mean=$(echo "${COMBATSEQ_TIMES[@]}" | tr ' ' '\n' | awk '{sum+=$1} END {print sum/NR}')
    combatseq_mean_s=$(echo "scale=2; $combatseq_mean / 1000" | bc)
    
    echo "combatseq activation: ${combatseq_mean_s}s (mean)"
    echo ""
else
    echo "⚠ No combatseq activation data collected"
    echo ""
fi

echo "Comparison with baseline (from task 12.1):"
echo "  Apptainer:     ~22ms"
echo "  uv (Python):   ~110ms"
echo "  nix-shell:     ~310,000ms (5 minutes)"
if [ "${#RESULT_TIMES[@]}" -gt 0 ] 2>/dev/null; then
    echo "  nix-shell opt: ~${result_mean_s}s (this test)"
else
    echo "  nix-shell opt: NOT TESTED (R environment not built)"
fi
echo ""

if [ "${#DEFAULT_TIMES[@]}" -gt 0 ] 2>/dev/null && [ "${#RESULT_TIMES[@]}" -gt 0 ] 2>/dev/null; then
    default_mean=$(echo "${DEFAULT_TIMES[@]}" | tr ' ' '\n' | awk '{sum+=$1} END {print sum/NR}')
    result_mean=$(echo "${RESULT_TIMES[@]}" | tr ' ' '\n' | awk '{sum+=$1} END {print sum/NR}')
    speedup=$(echo "scale=1; $default_mean / $result_mean" | bc)
    
    echo "Measured speedup: ${speedup}x faster"
    echo ""
fi

echo ""
echo "============================================"
echo "Next Steps"
echo "============================================"
if ! [ "${#RESULT_TIMES[@]}" -gt 0 ] 2>/dev/null; then
    echo ""
    echo "To complete this benchmark:"
    echo "  1. Build the R environment:"
    echo "     cd $BATCH_EFFECTS_DIR"
    echo "     ../../build_nix_env.sh"
    echo ""
    echo "  2. Re-run this benchmark:"
    echo "     bash environments/benchmarks/benchmark_nix_build_optimization.sh"
    echo ""
    echo "Note: The R environment build is currently blocked by Task 8.3"
    echo "      which is waiting for the nix-build to complete."
else
    echo ""
    echo "✓ Benchmark complete!"
    echo ""
    echo "Task 12.2 status: COMPLETE"
    echo "- Activation time measured: ${result_mean_s}s"
    echo "- Success criteria (<5s): $([ $(echo "$result_mean < 5000" | bc) -eq 1 ] && echo "PASS" || echo "FAIL")"
fi
echo ""

echo "Done!"
