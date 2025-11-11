#!/bin/bash
# SLURM Integration Validation Script
# Tests job submission and environment activation on compute nodes
#
# Usage:
#   bash environments/validate_slurm_integration.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(dirname "$SCRIPT_DIR")"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test results
TESTS_PASSED=0
TESTS_FAILED=0
JOB_IDS=()

print_header() {
    echo ""
    echo "======================================================================"
    echo "  $1"
    echo "======================================================================"
}

print_section() {
    echo ""
    echo "--- $1 ---"
}

print_success() {
    echo -e "${GREEN}✓${NC} $1"
    TESTS_PASSED=$((TESTS_PASSED + 1))
}

print_failure() {
    echo -e "${RED}✗${NC} $1"
    TESTS_FAILED=$((TESTS_FAILED + 1))
}

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

wait_for_job() {
    local job_id=$1
    local timeout=${2:-300}  # 5 minutes default
    local elapsed=0
    
    echo "  Waiting for job $job_id to complete (timeout: ${timeout}s)..."
    
    while [ $elapsed -lt $timeout ]; do
        # Check job status
        if ! squeue -j $job_id &>/dev/null; then
            # Job no longer in queue, check if it completed
            local state=$(sacct -j $job_id --format=State --noheader | head -1 | tr -d ' ')
            echo "  Job $job_id finished with state: $state"
            
            if [ "$state" = "COMPLETED" ]; then
                return 0
            else
                return 1
            fi
        fi
        
        sleep 5
        elapsed=$((elapsed + 5))
    done
    
    echo "  Job $job_id timed out after ${timeout}s"
    return 1
}

check_job_output() {
    local job_id=$1
    local expected_pattern=$2
    local output_file="slurm-${job_id}.out"
    
    if [ ! -f "$output_file" ]; then
        echo "  Output file not found: $output_file"
        return 1
    fi
    
    if grep -q "$expected_pattern" "$output_file"; then
        echo "  Found expected pattern in output: $expected_pattern"
        return 0
    else
        echo "  Expected pattern not found: $expected_pattern"
        echo "  Output file contents:"
        cat "$output_file" | head -20
        return 1
    fi
}

cleanup_job_files() {
    local job_id=$1
    rm -f "slurm-${job_id}.out" "slurm-${job_id}.err"
}

# ============================================================================
# Test 1: Python Script Submission
# ============================================================================
test_python_submission() {
    print_section "Test 1: Python Script SLURM Submission"
    
    # Submit Python validation script
    echo "  Submitting Python validation script..."
    local output=$(${SCRIPT_DIR}/run_with_env.sh --sbatch \
        --time 00:10:00 \
        --mem 4G \
        --cpus-per-task 2 \
        ${SCRIPT_DIR}/validate_python_imports.py 2>&1)
    
    # Extract job ID
    local job_id=$(echo "$output" | grep -oP 'Submitted batch job \K\d+' || echo "")
    
    if [ -z "$job_id" ]; then
        print_failure "Failed to submit Python job"
        echo "  Output: $output"
        return 1
    fi
    
    echo "  Job ID: $job_id"
    JOB_IDS+=($job_id)
    
    # Wait for job to complete
    if wait_for_job $job_id 600; then
        # Check output
        if check_job_output $job_id "ALL TESTS PASSED"; then
            print_success "Python script executed successfully on compute node"
            cleanup_job_files $job_id
            return 0
        else
            print_failure "Python script failed on compute node"
            return 1
        fi
    else
        print_failure "Python job did not complete successfully"
        return 1
    fi
}

# ============================================================================
# Test 2: R Script Submission (batch-effects environment)
# ============================================================================
test_r_batch_effects_submission() {
    print_section "Test 2: R Script SLURM Submission (batch-effects)"
    
    # Submit R validation script
    echo "  Submitting R validation script (batch-effects environment)..."
    local output=$(${SCRIPT_DIR}/run_with_env.sh --sbatch \
        --time 00:15:00 \
        --mem 8G \
        --cpus-per-task 2 \
        ${SCRIPT_DIR}/validate_r_packages.R 2>&1)
    
    # Extract job ID
    local job_id=$(echo "$output" | grep -oP 'Submitted batch job \K\d+' || echo "")
    
    if [ -z "$job_id" ]; then
        print_failure "Failed to submit R job (batch-effects)"
        echo "  Output: $output"
        return 1
    fi
    
    echo "  Job ID: $job_id"
    JOB_IDS+=($job_id)
    
    # Wait for job to complete
    if wait_for_job $job_id 900; then
        # Check output
        if check_job_output $job_id "ALL TESTS PASSED"; then
            print_success "R script (batch-effects) executed successfully on compute node"
            cleanup_job_files $job_id
            return 0
        else
            print_failure "R script (batch-effects) failed on compute node"
            return 1
        fi
    else
        print_failure "R job (batch-effects) did not complete successfully"
        return 1
    fi
}

# ============================================================================
# Test 3: R Script Submission (combatseq environment)
# ============================================================================
test_r_combatseq_submission() {
    print_section "Test 3: R Script SLURM Submission (combatseq)"
    
    # Check if combatseq environment exists
    if [ ! -d "${ANALYSIS_DIR}/environments/r/combatseq" ]; then
        print_warning "ComBat-seq environment not found, skipping test"
        return 0
    fi
    
    # Submit R validation script with combatseq environment
    echo "  Submitting R validation script (combatseq environment)..."
    local output=$(${SCRIPT_DIR}/run_with_env.sh --r-env combatseq --sbatch \
        --time 00:15:00 \
        --mem 8G \
        --cpus-per-task 2 \
        ${SCRIPT_DIR}/validate_r_packages.R 2>&1)
    
    # Extract job ID
    local job_id=$(echo "$output" | grep -oP 'Submitted batch job \K\d+' || echo "")
    
    if [ -z "$job_id" ]; then
        print_failure "Failed to submit R job (combatseq)"
        echo "  Output: $output"
        return 1
    fi
    
    echo "  Job ID: $job_id"
    JOB_IDS+=($job_id)
    
    # Wait for job to complete
    if wait_for_job $job_id 900; then
        # Check output
        if check_job_output $job_id "ALL TESTS PASSED"; then
            print_success "R script (combatseq) executed successfully on compute node"
            cleanup_job_files $job_id
            return 0
        else
            print_failure "R script (combatseq) failed on compute node"
            return 1
        fi
    else
        print_failure "R job (combatseq) did not complete successfully"
        return 1
    fi
}

# ============================================================================
# Test 4: Mixed Environment Script
# ============================================================================
test_mixed_environment() {
    print_section "Test 4: Mixed Environment (Python + R)"
    
    # Create a simple test script that uses both Python and R
    local test_script=$(mktemp /tmp/test_mixed_env.XXXXXX.R)
    cat > "$test_script" << 'EOF'
#!/usr/bin/env Rscript
# Test script that uses both R and Python (via reticulate)

cat("Testing mixed environment (R + Python)\n")

# Test R functionality
cat("Testing R...\n")
library(dplyr)
df <- data.frame(x = 1:10, y = rnorm(10))
result <- df %>% summarise(mean_y = mean(y))
cat(sprintf("  R: mean(y) = %.4f\n", result$mean_y))

# Test Python via reticulate (if available)
if (requireNamespace("reticulate", quietly = TRUE)) {
  cat("Testing Python via reticulate...\n")
  library(reticulate)
  
  # Check if Python is available
  py_available <- tryCatch({
    py_config()
    TRUE
  }, error = function(e) {
    cat(sprintf("  Warning: Python not available: %s\n", e$message))
    FALSE
  })
  
  if (py_available) {
    # Test Python import
    np <- import("numpy", convert = FALSE)
    py_arr <- np$array(c(1, 2, 3, 4, 5))
    py_mean <- np$mean(py_arr)
    cat(sprintf("  Python: numpy.mean([1,2,3,4,5]) = %.4f\n", py_mean))
    cat("  ✓ Mixed environment test passed\n")
  }
} else {
  cat("  Warning: reticulate not available, skipping Python test\n")
}

cat("Mixed environment test completed\n")
EOF
    
    chmod +x "$test_script"
    
    # Submit with full environment
    echo "  Submitting mixed environment test..."
    local output=$(${SCRIPT_DIR}/run_with_env.sh --full-env --sbatch \
        --time 00:10:00 \
        --mem 4G \
        --cpus-per-task 2 \
        "$test_script" 2>&1)
    
    # Extract job ID
    local job_id=$(echo "$output" | grep -oP 'Submitted batch job \K\d+' || echo "")
    
    if [ -z "$job_id" ]; then
        print_failure "Failed to submit mixed environment job"
        echo "  Output: $output"
        rm -f "$test_script"
        return 1
    fi
    
    echo "  Job ID: $job_id"
    JOB_IDS+=($job_id)
    
    # Wait for job to complete
    if wait_for_job $job_id 600; then
        # Check output
        if check_job_output $job_id "Mixed environment test completed"; then
            print_success "Mixed environment script executed successfully"
            cleanup_job_files $job_id
            rm -f "$test_script"
            return 0
        else
            print_failure "Mixed environment script failed"
            rm -f "$test_script"
            return 1
        fi
    else
        print_failure "Mixed environment job did not complete successfully"
        rm -f "$test_script"
        return 1
    fi
}

# ============================================================================
# Test 5: Array Job Submission
# ============================================================================
test_array_job() {
    print_section "Test 5: Array Job Submission"
    
    # Create a simple array job test script
    local test_script=$(mktemp /tmp/test_array_job.XXXXXX.py)
    cat > "$test_script" << 'EOF'
#!/usr/bin/env python3
import os
import sys

task_id = os.environ.get('SLURM_ARRAY_TASK_ID', 'unknown')
job_id = os.environ.get('SLURM_ARRAY_JOB_ID', 'unknown')

print(f"Array job test: Job {job_id}, Task {task_id}")
print(f"Python version: {sys.version.split()[0]}")

# Simple computation
import numpy as np
result = np.sum(np.arange(int(task_id) * 100))
print(f"Task {task_id} result: {result}")
print("Array task completed successfully")
EOF
    
    chmod +x "$test_script"
    
    # Submit array job
    echo "  Submitting array job (5 tasks)..."
    local output=$(${SCRIPT_DIR}/run_with_env.sh --sbatch \
        --array 1-5 \
        --time 00:05:00 \
        --mem 2G \
        --cpus-per-task 1 \
        "$test_script" 2>&1)
    
    # Extract job ID
    local job_id=$(echo "$output" | grep -oP 'Submitted batch job \K\d+' || echo "")
    
    if [ -z "$job_id" ]; then
        print_failure "Failed to submit array job"
        echo "  Output: $output"
        rm -f "$test_script"
        return 1
    fi
    
    echo "  Job ID: $job_id"
    JOB_IDS+=($job_id)
    
    # Wait for all array tasks to complete
    if wait_for_job $job_id 600; then
        # Check that all 5 tasks completed
        local completed_tasks=$(sacct -j $job_id --format=JobID,State --noheader | grep -c "COMPLETED" || echo "0")
        
        if [ "$completed_tasks" -ge 5 ]; then
            print_success "Array job completed successfully (${completed_tasks} tasks)"
            
            # Clean up array job output files
            for i in {1..5}; do
                rm -f "slurm-${job_id}_${i}.out" "slurm-${job_id}_${i}.err"
            done
            rm -f "$test_script"
            return 0
        else
            print_failure "Array job incomplete (only ${completed_tasks}/5 tasks completed)"
            rm -f "$test_script"
            return 1
        fi
    else
        print_failure "Array job did not complete successfully"
        rm -f "$test_script"
        return 1
    fi
}

# ============================================================================
# Test 6: Environment Activation Check
# ============================================================================
test_environment_activation() {
    print_section "Test 6: Environment Activation on Compute Node"
    
    # Create a script that checks environment variables
    local test_script=$(mktemp /tmp/test_env_activation.XXXXXX.sh)
    cat > "$test_script" << 'EOF'
#!/bin/bash
set -euo pipefail

echo "=== Environment Activation Test ==="
echo ""
echo "Hostname: $(hostname)"
echo "Node type: $(hostname | grep -q 'compute\|node\|cn' && echo 'Compute Node' || echo 'Login Node')"
echo ""

# Check Python environment
echo "Python Environment:"
echo "  PYTHON_ENV: ${PYTHON_ENV:-not set}"
echo "  Python executable: $(which python)"
echo "  Python version: $(python --version)"
echo ""

# Check R environment
echo "R Environment:"
echo "  R_ENV_DIR: ${R_ENV_DIR:-not set}"
echo "  R executable: $(which R || echo 'not found')"
if command -v R &>/dev/null; then
    echo "  R version: $(R --version | head -1)"
fi
echo ""

# Test Python import
echo "Testing Python import..."
python -c "import numpy; print(f'  NumPy version: {numpy.__version__}')"
echo ""

# Test R library
if command -v R &>/dev/null; then
    echo "Testing R library..."
    Rscript -e "library(dplyr); cat('  dplyr loaded successfully\n')"
    echo ""
fi

echo "✓ Environment activation test completed"
EOF
    
    chmod +x "$test_script"
    
    # Submit test script
    echo "  Submitting environment activation test..."
    local output=$(${SCRIPT_DIR}/run_with_env.sh --sbatch \
        --time 00:05:00 \
        --mem 2G \
        --cpus-per-task 1 \
        "$test_script" 2>&1)
    
    # Extract job ID
    local job_id=$(echo "$output" | grep -oP 'Submitted batch job \K\d+' || echo "")
    
    if [ -z "$job_id" ]; then
        print_failure "Failed to submit environment activation test"
        echo "  Output: $output"
        rm -f "$test_script"
        return 1
    fi
    
    echo "  Job ID: $job_id"
    JOB_IDS+=($job_id)
    
    # Wait for job to complete
    if wait_for_job $job_id 300; then
        # Check output
        if check_job_output $job_id "Environment activation test completed"; then
            print_success "Environment activated successfully on compute node"
            
            # Show relevant output
            echo "  Environment details:"
            grep -A 2 "Python Environment:" "slurm-${job_id}.out" || true
            grep -A 2 "R Environment:" "slurm-${job_id}.out" || true
            
            cleanup_job_files $job_id
            rm -f "$test_script"
            return 0
        else
            print_failure "Environment activation test failed"
            rm -f "$test_script"
            return 1
        fi
    else
        print_failure "Environment activation test did not complete successfully"
        rm -f "$test_script"
        return 1
    fi
}

# ============================================================================
# Main Function
# ============================================================================
main() {
    print_header "SLURM Integration Validation"
    echo "  Started: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "  Working directory: $(pwd)"
    echo ""
    
    # Check prerequisites
    if ! command -v sbatch &>/dev/null; then
        echo "Error: sbatch command not found. Are you on an HPC cluster?"
        exit 1
    fi
    
    if [ ! -f "${SCRIPT_DIR}/run_with_env.sh" ]; then
        echo "Error: run_with_env.sh not found at ${SCRIPT_DIR}/run_with_env.sh"
        exit 1
    fi
    
    # Run tests
    test_python_submission || true
    test_r_batch_effects_submission || true
    test_r_combatseq_submission || true
    test_mixed_environment || true
    test_array_job || true
    test_environment_activation || true
    
    # Summary
    print_header "Test Summary"
    
    local total_tests=$((TESTS_PASSED + TESTS_FAILED))
    echo ""
    echo "  Total tests:   $total_tests"
    echo "  Passed:        $TESTS_PASSED"
    echo "  Failed:        $TESTS_FAILED"
    echo ""
    
    if [ ${#JOB_IDS[@]} -gt 0 ]; then
        echo "  Submitted job IDs: ${JOB_IDS[*]}"
        echo ""
    fi
    
    if [ $TESTS_FAILED -eq 0 ]; then
        echo -e "${GREEN}✓ ALL SLURM INTEGRATION TESTS PASSED${NC}"
        echo "======================================================================"
        return 0
    else
        echo -e "${RED}✗ SOME SLURM INTEGRATION TESTS FAILED${NC}"
        echo "======================================================================"
        return 1
    fi
}

# Run main function
main
exit $?
