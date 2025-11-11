#!/bin/bash
# benchmark_array_jobs.sh
# Benchmark SLURM array jobs for Apptainer vs uv/rix
#
# Tests:
# - Small array job (10 jobs) with uv/rix
# - Medium array job (100 jobs) with uv/rix
# - Compare with equivalent Apptainer array jobs
# - Test with both Python and R scripts

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
RESULTS_DIR="${SCRIPT_DIR}/array_job_results"

mkdir -p "$RESULTS_DIR"

echo "=== Benchmark: SLURM Array Jobs ==="
echo "Results directory: $RESULTS_DIR"
echo ""

# ============================================================================
# Create Test Scripts
# ============================================================================

# Python test script
cat > "${RESULTS_DIR}/test_python.py" << 'EOF'
#!/usr/bin/env python
"""Simple Python test script for array job benchmarking."""
import sys
import time
import numpy as np

# Get array task ID from environment or argument
task_id = int(sys.argv[1]) if len(sys.argv) > 1 else 0

# Simulate some work
start = time.time()
data = np.random.rand(1000, 1000)
result = np.sum(data)
elapsed = time.time() - start

print(f"Task {task_id}: Computed sum = {result:.2f} in {elapsed:.3f}s")
EOF

chmod +x "${RESULTS_DIR}/test_python.py"

# R test script
cat > "${RESULTS_DIR}/test_r.R" << 'EOF'
#!/usr/bin/env Rscript
# Simple R test script for array job benchmarking

# Get array task ID from environment or argument
args <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) > 0) as.integer(args[1]) else 0

# Simulate some work
start_time <- Sys.time()
data <- matrix(rnorm(1000 * 1000), nrow = 1000)
result <- sum(data)
elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat(sprintf("Task %d: Computed sum = %.2f in %.3fs\n", task_id, result, elapsed))
EOF

chmod +x "${RESULTS_DIR}/test_r.R"

echo "Created test scripts:"
echo "  - ${RESULTS_DIR}/test_python.py"
echo "  - ${RESULTS_DIR}/test_r.R"
echo ""

# ============================================================================
# Helper Functions
# ============================================================================

# Wait for job to complete and collect statistics
wait_for_job() {
    local job_id="$1"
    local job_name="$2"
    
    echo "Waiting for job $job_id ($job_name) to complete..."
    
    # Wait for job to finish
    while squeue -j "$job_id" 2>/dev/null | grep -q "$job_id"; do
        sleep 5
    done
    
    echo "Job $job_id completed. Collecting statistics..."
    
    # Get job statistics
    sacct -j "$job_id" --format=JobID,State,ExitCode,Elapsed,MaxRSS,CPUTime --parsable2 > "${RESULTS_DIR}/${job_name}_stats.txt"
    
    # Count successful and failed tasks
    local total=$(sacct -j "$job_id" --format=JobID,State --parsable2 | grep -c "_" || true)
    local completed=$(sacct -j "$job_id" --format=JobID,State --parsable2 | grep "_" | grep -c "COMPLETED" || true)
    local failed=$(sacct -j "$job_id" --format=JobID,State --parsable2 | grep "_" | grep -c "FAILED" || true)
    
    # Calculate total wall time (from first task start to last task end)
    local start_time=$(sacct -j "$job_id" --format=Start --parsable2 | tail -n +2 | head -n 1)
    local end_time=$(sacct -j "$job_id" --format=End --parsable2 | tail -n +2 | sort | tail -n 1)
    
    echo ""
    echo "Job $job_id ($job_name) Summary:"
    echo "  Total tasks:     $total"
    echo "  Completed:       $completed"
    echo "  Failed:          $failed"
    echo "  Start time:      $start_time"
    echo "  End time:        $end_time"
    echo "  Statistics saved to: ${RESULTS_DIR}/${job_name}_stats.txt"
    echo ""
}

# ============================================================================
# Benchmark 1: Small Array Job (10 tasks) - Python
# ============================================================================
echo "=== Benchmark 1: Small Array Job (10 tasks) - Python ==="
echo ""

# Apptainer version
echo "Submitting Apptainer array job (10 tasks, Python)..."
APPTAINER_SMALL_PY_JOB=$(sbatch --parsable \
    --job-name=bench_app_small_py \
    --array=1-10 \
    --time=00:05:00 \
    --mem=1G \
    --cpus-per-task=1 \
    --output="${RESULTS_DIR}/app_small_py_%a.out" \
    --wrap="sg grp_batch_effects -c 'apptainer exec --contain ~/groups/grp_batch_effects/remove-batch-effects.sif python ${RESULTS_DIR}/test_python.py \$SLURM_ARRAY_TASK_ID'")

echo "Submitted job: $APPTAINER_SMALL_PY_JOB"

# uv version
echo "Submitting uv array job (10 tasks, Python)..."
source "${ANALYSIS_DIR}/environments/init_env.sh" 2>/dev/null

UV_SMALL_PY_JOB=$(sbatch --parsable \
    --job-name=bench_uv_small_py \
    --array=1-10 \
    --time=00:05:00 \
    --mem=1G \
    --cpus-per-task=1 \
    --output="${RESULTS_DIR}/uv_small_py_%a.out" \
    --wrap="source ${ANALYSIS_DIR}/environments/init_env.sh && source ${PYTHON_ENV}/bin/activate && python ${RESULTS_DIR}/test_python.py \$SLURM_ARRAY_TASK_ID")

echo "Submitted job: $UV_SMALL_PY_JOB"
echo ""

# ============================================================================
# Benchmark 2: Small Array Job (10 tasks) - R
# ============================================================================
echo "=== Benchmark 2: Small Array Job (10 tasks) - R ==="
echo ""

NIX_ROOT="/grphome/grp_batch_effects/nix"
R_ENV_DIR="${ANALYSIS_DIR}/environments/r/batch-effects"
NIX_CACHE_OPTS="--option substituters 'https://cache.nixos.org https://rstats-on-nix.cachix.org' --option trusted-public-keys 'cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= rstats-on-nix.cachix.org-1:9cJb6nqYZgKqgH5XJQN8FPkXqKlGqKqJqKqKqKqKqKo='"

# Check if R environment exists
if [ ! -d "$R_ENV_DIR" ] || [ ! -f "$R_ENV_DIR/default.nix" ]; then
    echo "WARNING: R environment not found at $R_ENV_DIR"
    echo "Skipping R array job benchmarks. Run Phase 1 (authoring) first."
    echo ""
    R_BENCHMARKS_SKIPPED=true
else
    R_BENCHMARKS_SKIPPED=false
    
    # Apptainer version
    echo "Submitting Apptainer array job (10 tasks, R)..."
    APPTAINER_SMALL_R_JOB=$(sbatch --parsable \
        --job-name=bench_app_small_r \
        --array=1-10 \
        --time=00:05:00 \
        --mem=1G \
        --cpus-per-task=1 \
        --output="${RESULTS_DIR}/app_small_r_%a.out" \
        --wrap="sg grp_batch_effects -c 'apptainer exec --contain ~/groups/grp_batch_effects/remove-batch-effects.sif Rscript ${RESULTS_DIR}/test_r.R \$SLURM_ARRAY_TASK_ID'")
    
    echo "Submitted job: $APPTAINER_SMALL_R_JOB"
    
    # rix version
    echo "Submitting rix array job (10 tasks, R)..."
    RIX_SMALL_R_JOB=$(sbatch --parsable \
        --job-name=bench_rix_small_r \
        --array=1-10 \
        --time=00:05:00 \
        --mem=1G \
        --cpus-per-task=1 \
        --output="${RESULTS_DIR}/rix_small_r_%a.out" \
        --wrap="$NIX_ROOT/nix-user-chroot $NIX_ROOT bash -c 'source ~/.nix-profile/etc/profile.d/nix.sh && cd $R_ENV_DIR && nix-shell $NIX_CACHE_OPTS --run \"Rscript ${RESULTS_DIR}/test_r.R \$SLURM_ARRAY_TASK_ID\"'")
    
    echo "Submitted job: $RIX_SMALL_R_JOB"
    echo ""
fi

# ============================================================================
# Benchmark 3: Medium Array Job (100 tasks) - Python
# ============================================================================
echo "=== Benchmark 3: Medium Array Job (100 tasks) - Python ==="
echo ""

# Apptainer version
echo "Submitting Apptainer array job (100 tasks, Python)..."
APPTAINER_MEDIUM_PY_JOB=$(sbatch --parsable \
    --job-name=bench_app_medium_py \
    --array=1-100 \
    --time=00:05:00 \
    --mem=1G \
    --cpus-per-task=1 \
    --output="${RESULTS_DIR}/app_medium_py_%a.out" \
    --wrap="sg grp_batch_effects -c 'apptainer exec --contain ~/groups/grp_batch_effects/remove-batch-effects.sif python ${RESULTS_DIR}/test_python.py \$SLURM_ARRAY_TASK_ID'")

echo "Submitted job: $APPTAINER_MEDIUM_PY_JOB"

# uv version
echo "Submitting uv array job (100 tasks, Python)..."
UV_MEDIUM_PY_JOB=$(sbatch --parsable \
    --job-name=bench_uv_medium_py \
    --array=1-100 \
    --time=00:05:00 \
    --mem=1G \
    --cpus-per-task=1 \
    --output="${RESULTS_DIR}/uv_medium_py_%a.out" \
    --wrap="source ${ANALYSIS_DIR}/environments/init_env.sh && source ${PYTHON_ENV}/bin/activate && python ${RESULTS_DIR}/test_python.py \$SLURM_ARRAY_TASK_ID")

echo "Submitted job: $UV_MEDIUM_PY_JOB"
echo ""

# ============================================================================
# Benchmark 4: Medium Array Job (100 tasks) - R
# ============================================================================
if [ "$R_BENCHMARKS_SKIPPED" = false ]; then
    echo "=== Benchmark 4: Medium Array Job (100 tasks) - R ==="
    echo ""
    
    # Apptainer version
    echo "Submitting Apptainer array job (100 tasks, R)..."
    APPTAINER_MEDIUM_R_JOB=$(sbatch --parsable \
        --job-name=bench_app_medium_r \
        --array=1-100 \
        --time=00:05:00 \
        --mem=1G \
        --cpus-per-task=1 \
        --output="${RESULTS_DIR}/app_medium_r_%a.out" \
        --wrap="sg grp_batch_effects -c 'apptainer exec --contain ~/groups/grp_batch_effects/remove-batch-effects.sif Rscript ${RESULTS_DIR}/test_r.R \$SLURM_ARRAY_TASK_ID'")
    
    echo "Submitted job: $APPTAINER_MEDIUM_R_JOB"
    
    # rix version
    echo "Submitting rix array job (100 tasks, R)..."
    RIX_MEDIUM_R_JOB=$(sbatch --parsable \
        --job-name=bench_rix_medium_r \
        --array=1-100 \
        --time=00:05:00 \
        --mem=1G \
        --cpus-per-task=1 \
        --output="${RESULTS_DIR}/rix_medium_r_%a.out" \
        --wrap="$NIX_ROOT/nix-user-chroot $NIX_ROOT bash -c 'source ~/.nix-profile/etc/profile.d/nix.sh && cd $R_ENV_DIR && nix-shell $NIX_CACHE_OPTS --run \"Rscript ${RESULTS_DIR}/test_r.R \$SLURM_ARRAY_TASK_ID\"'")
    
    echo "Submitted job: $RIX_MEDIUM_R_JOB"
    echo ""
fi

# ============================================================================
# Wait for Jobs and Collect Statistics
# ============================================================================
echo "=== Waiting for Jobs to Complete ==="
echo ""
echo "This may take several minutes. You can monitor jobs with:"
echo "  squeue -u \$USER"
echo ""

# Wait for all jobs
wait_for_job "$APPTAINER_SMALL_PY_JOB" "apptainer_small_py"
wait_for_job "$UV_SMALL_PY_JOB" "uv_small_py"

if [ "$R_BENCHMARKS_SKIPPED" = false ]; then
    wait_for_job "$APPTAINER_SMALL_R_JOB" "apptainer_small_r"
    wait_for_job "$RIX_SMALL_R_JOB" "rix_small_r"
fi

wait_for_job "$APPTAINER_MEDIUM_PY_JOB" "apptainer_medium_py"
wait_for_job "$UV_MEDIUM_PY_JOB" "uv_medium_py"

if [ "$R_BENCHMARKS_SKIPPED" = false ]; then
    wait_for_job "$APPTAINER_MEDIUM_R_JOB" "apptainer_medium_r"
    wait_for_job "$RIX_MEDIUM_R_JOB" "rix_medium_r"
fi

# ============================================================================
# Generate Summary Report
# ============================================================================
echo "=== Generating Summary Report ==="
echo ""

SUMMARY_FILE="${RESULTS_DIR}/summary.txt"

cat > "$SUMMARY_FILE" << EOF
SLURM Array Job Benchmark Summary
==================================

Test Date: $(date)

Small Array Jobs (10 tasks):
-----------------------------
EOF

# Python small
if [ -f "${RESULTS_DIR}/apptainer_small_py_stats.txt" ]; then
    echo "" >> "$SUMMARY_FILE"
    echo "Python - Apptainer (10 tasks):" >> "$SUMMARY_FILE"
    tail -n +2 "${RESULTS_DIR}/apptainer_small_py_stats.txt" | head -n 5 >> "$SUMMARY_FILE"
fi

if [ -f "${RESULTS_DIR}/uv_small_py_stats.txt" ]; then
    echo "" >> "$SUMMARY_FILE"
    echo "Python - uv (10 tasks):" >> "$SUMMARY_FILE"
    tail -n +2 "${RESULTS_DIR}/uv_small_py_stats.txt" | head -n 5 >> "$SUMMARY_FILE"
fi

# R small
if [ "$R_BENCHMARKS_SKIPPED" = false ]; then
    if [ -f "${RESULTS_DIR}/apptainer_small_r_stats.txt" ]; then
        echo "" >> "$SUMMARY_FILE"
        echo "R - Apptainer (10 tasks):" >> "$SUMMARY_FILE"
        tail -n +2 "${RESULTS_DIR}/apptainer_small_r_stats.txt" | head -n 5 >> "$SUMMARY_FILE"
    fi
    
    if [ -f "${RESULTS_DIR}/rix_small_r_stats.txt" ]; then
        echo "" >> "$SUMMARY_FILE"
        echo "R - rix (10 tasks):" >> "$SUMMARY_FILE"
        tail -n +2 "${RESULTS_DIR}/rix_small_r_stats.txt" | head -n 5 >> "$SUMMARY_FILE"
    fi
fi

cat >> "$SUMMARY_FILE" << EOF

Medium Array Jobs (100 tasks):
-------------------------------
EOF

# Python medium
if [ -f "${RESULTS_DIR}/apptainer_medium_py_stats.txt" ]; then
    echo "" >> "$SUMMARY_FILE"
    echo "Python - Apptainer (100 tasks):" >> "$SUMMARY_FILE"
    tail -n +2 "${RESULTS_DIR}/apptainer_medium_py_stats.txt" | head -n 5 >> "$SUMMARY_FILE"
fi

if [ -f "${RESULTS_DIR}/uv_medium_py_stats.txt" ]; then
    echo "" >> "$SUMMARY_FILE"
    echo "Python - uv (100 tasks):" >> "$SUMMARY_FILE"
    tail -n +2 "${RESULTS_DIR}/uv_medium_py_stats.txt" | head -n 5 >> "$SUMMARY_FILE"
fi

# R medium
if [ "$R_BENCHMARKS_SKIPPED" = false ]; then
    if [ -f "${RESULTS_DIR}/apptainer_medium_r_stats.txt" ]; then
        echo "" >> "$SUMMARY_FILE"
        echo "R - Apptainer (100 tasks):" >> "$SUMMARY_FILE"
        tail -n +2 "${RESULTS_DIR}/apptainer_medium_r_stats.txt" | head -n 5 >> "$SUMMARY_FILE"
    fi
    
    if [ -f "${RESULTS_DIR}/rix_medium_r_stats.txt" ]; then
        echo "" >> "$SUMMARY_FILE"
        echo "R - rix (100 tasks):" >> "$SUMMARY_FILE"
        tail -n +2 "${RESULTS_DIR}/rix_medium_r_stats.txt" | head -n 5 >> "$SUMMARY_FILE"
    fi
fi

echo "" >> "$SUMMARY_FILE"
echo "All results saved to: $RESULTS_DIR" >> "$SUMMARY_FILE"

# Display summary
cat "$SUMMARY_FILE"

echo ""
echo "=== Benchmark Complete ==="
echo ""
echo "Results directory: $RESULTS_DIR"
echo "Summary report: $SUMMARY_FILE"
echo ""
echo "Individual job outputs:"
ls -lh "${RESULTS_DIR}"/*.out 2>/dev/null || echo "  (No output files found)"
echo ""
echo "Job statistics:"
ls -lh "${RESULTS_DIR}"/*_stats.txt 2>/dev/null || echo "  (No stats files found)"
echo ""
