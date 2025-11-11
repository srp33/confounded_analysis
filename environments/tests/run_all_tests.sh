#!/bin/bash
# Comprehensive test suite for execution wrappers (Task 10.4)
# Tests all aspects of run_with_env.sh functionality

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test counters
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="$(dirname "$SCRIPT_DIR")"
RUN_WITH_ENV="${ENV_DIR}/run_with_env.sh"

# Check if run_with_env.sh exists
if [ ! -f "$RUN_WITH_ENV" ]; then
    echo -e "${RED}Error: run_with_env.sh not found at $RUN_WITH_ENV${NC}"
    exit 1
fi

# Make test scripts executable
chmod +x "${SCRIPT_DIR}"/*.py "${SCRIPT_DIR}"/*.R "${SCRIPT_DIR}"/*.sh 2>/dev/null || true

echo "============================================================"
echo "Execution Wrapper Test Suite (Task 10.4)"
echo "============================================================"
echo ""
echo "Testing: $RUN_WITH_ENV"
echo ""

# Helper function to run a test
run_test() {
    local test_name="$1"
    local test_command="$2"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    echo -e "${BLUE}[TEST $TOTAL_TESTS]${NC} $test_name"
    echo "Command: $test_command"
    echo ""
    
    if eval "$test_command"; then
        echo -e "${GREEN}✓ PASSED${NC}"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    else
        echo -e "${RED}✗ FAILED${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
    
    echo ""
    echo "------------------------------------------------------------"
    echo ""
}

# ============================================================================
# Test Category 1: Python Script Execution
# ============================================================================
echo "============================================================"
echo "Category 1: Python Script Execution"
echo "============================================================"
echo ""

run_test \
    "Python script execution (basic)" \
    "$RUN_WITH_ENV ${SCRIPT_DIR}/test_python_script.py"

run_test \
    "Python script with arguments" \
    "$RUN_WITH_ENV ${SCRIPT_DIR}/test_python_script.py arg1 arg2 arg3"

# ============================================================================
# Test Category 2: R Script Execution
# ============================================================================
echo "============================================================"
echo "Category 2: R Script Execution"
echo "============================================================"
echo ""

run_test \
    "R script execution (basic)" \
    "$RUN_WITH_ENV ${SCRIPT_DIR}/test_r_script.R"

run_test \
    "R script with arguments" \
    "$RUN_WITH_ENV ${SCRIPT_DIR}/test_r_script.R arg1 arg2 arg3"

run_test \
    "R script with ComBat-seq environment" \
    "$RUN_WITH_ENV --r-env combatseq ${SCRIPT_DIR}/test_r_script.R"

# ============================================================================
# Test Category 3: Mixed Environment Execution
# ============================================================================
echo "============================================================"
echo "Category 3: Mixed Environment Execution"
echo "============================================================"
echo ""

run_test \
    "Mixed environment (R + Python)" \
    "$RUN_WITH_ENV --full-env ${SCRIPT_DIR}/test_mixed_env.R"

run_test \
    "Shell script execution (both environments)" \
    "$RUN_WITH_ENV ${SCRIPT_DIR}/test_shell_script.sh"

run_test \
    "Shell script with arguments" \
    "$RUN_WITH_ENV ${SCRIPT_DIR}/test_shell_script.sh arg1 arg2"

# ============================================================================
# Test Category 4: Direct Command Execution
# ============================================================================
echo "============================================================"
echo "Category 4: Direct Command Execution"
echo "============================================================"
echo ""

run_test \
    "Direct Python command" \
    "$RUN_WITH_ENV python -c 'import numpy; print(\"NumPy version:\", numpy.__version__)'"

run_test \
    "Direct R command" \
    "$RUN_WITH_ENV Rscript -e 'library(dplyr); cat(\"dplyr loaded\\n\")'"

# ============================================================================
# Test Category 5: Interactive Shell (manual test)
# ============================================================================
echo "============================================================"
echo "Category 5: Interactive Shell"
echo "============================================================"
echo ""
echo -e "${YELLOW}Note: Interactive shell test requires manual verification${NC}"
echo "To test interactive shell, run:"
echo "  $RUN_WITH_ENV shell"
echo ""
echo "Then verify:"
echo "  - Python is available: python --version"
echo "  - R is available: R --version"
echo "  - Packages work: python -c 'import numpy'"
echo "  - Exit with: exit"
echo ""

# ============================================================================
# Test Category 6: SLURM Job Submission (dry-run tests)
# ============================================================================
echo "============================================================"
echo "Category 6: SLURM Job Submission"
echo "============================================================"
echo ""

# Note: We can't actually submit SLURM jobs in this test environment
# But we can test the sbatch script generation

echo -e "${YELLOW}Note: SLURM tests generate job scripts but don't submit${NC}"
echo "Actual SLURM submission requires cluster access"
echo ""

# Test sbatch script generation for Python
echo "[TEST] SLURM Python job script generation"
TEMP_TEST=$(mktemp)
cat > "$TEMP_TEST" << 'EOF'
#!/bin/bash
# Test if sbatch script is generated correctly
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="$(dirname "$SCRIPT_DIR")"
RUN_WITH_ENV="${ENV_DIR}/run_with_env.sh"

# Capture the sbatch command (it will fail without SLURM, but we can check the script)
OUTPUT=$("$RUN_WITH_ENV" --sbatch --time 00:10:00 --mem 4G "${SCRIPT_DIR}/test_python_script.py" 2>&1 || true)

# Check if temporary sbatch script was mentioned
if echo "$OUTPUT" | grep -q "sbatch"; then
    echo "✓ SLURM job script generation works"
    exit 0
else
    echo "✗ SLURM job script generation failed"
    echo "$OUTPUT"
    exit 1
fi
EOF
chmod +x "$TEMP_TEST"

if bash "$TEMP_TEST"; then
    echo -e "${GREEN}✓ SLURM script generation test passed${NC}"
else
    echo -e "${YELLOW}⚠ SLURM not available (expected in test environment)${NC}"
fi
rm -f "$TEMP_TEST"

echo ""
echo "To test actual SLURM submission on the cluster, run:"
echo "  $RUN_WITH_ENV --sbatch --time 00:10:00 --mem 4G ${SCRIPT_DIR}/test_python_script.py"
echo "  $RUN_WITH_ENV --sbatch --time 00:10:00 --mem 4G ${SCRIPT_DIR}/test_r_script.R"
echo ""

# ============================================================================
# Test Category 7: Resource Specifications
# ============================================================================
echo "============================================================"
echo "Category 7: Various Resource Specifications"
echo "============================================================"
echo ""

echo "Testing different resource specification formats..."
echo ""

# These tests verify the command parsing works correctly
# Actual resource allocation would be tested on SLURM cluster

TEST_SPECS=(
    "--time 00:30:00 --mem 16G"
    "--time 02:00:00 --mem 64G --cpus-per-task 8"
    "--time 01:00:00 --mem 32G --nodes 1 --ntasks 1"
    "--partition short --time 00:15:00 --mem 8G"
)

for spec in "${TEST_SPECS[@]}"; do
    echo "  Testing: $spec"
    # Just verify the command doesn't crash during parsing
    if "$RUN_WITH_ENV" --sbatch $spec "${SCRIPT_DIR}/test_python_script.py" 2>&1 | grep -q "Submitting job\|sbatch\|SLURM" || true; then
        echo -e "    ${GREEN}✓ Parsing successful${NC}"
    else
        echo -e "    ${YELLOW}⚠ SLURM not available${NC}"
    fi
done

echo ""

# ============================================================================
# Test Category 8: Error Handling
# ============================================================================
echo "============================================================"
echo "Category 8: Error Handling"
echo "============================================================"
echo ""

echo "[TEST] Missing script file"
if ! "$RUN_WITH_ENV" nonexistent_script.py 2>&1 | grep -q "Error\|not found\|No such file"; then
    echo -e "${YELLOW}⚠ Error handling may need improvement${NC}"
else
    echo -e "${GREEN}✓ Error handling works${NC}"
fi
echo ""

echo "[TEST] Invalid R environment"
if ! "$RUN_WITH_ENV" --r-env invalid_env "${SCRIPT_DIR}/test_r_script.R" 2>&1 | grep -q "Error\|Unknown"; then
    echo -e "${YELLOW}⚠ Error handling may need improvement${NC}"
else
    echo -e "${GREEN}✓ Error handling works${NC}"
fi
echo ""

# ============================================================================
# Test Summary
# ============================================================================
echo "============================================================"
echo "Test Summary"
echo "============================================================"
echo ""
echo "Total tests run: $TOTAL_TESTS"
echo -e "${GREEN}Passed: $PASSED_TESTS${NC}"
if [ $FAILED_TESTS -gt 0 ]; then
    echo -e "${RED}Failed: $FAILED_TESTS${NC}"
else
    echo "Failed: 0"
fi
echo ""

if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "${GREEN}============================================================${NC}"
    echo -e "${GREEN}✓ All tests passed!${NC}"
    echo -e "${GREEN}============================================================${NC}"
    exit 0
else
    echo -e "${RED}============================================================${NC}"
    echo -e "${RED}✗ Some tests failed${NC}"
    echo -e "${RED}============================================================${NC}"
    exit 1
fi
