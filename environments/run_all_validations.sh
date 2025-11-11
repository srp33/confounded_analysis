#!/bin/bash
# Master Validation Script
# Runs all validation tests for Task 11 (Library import validation)
#
# Usage:
#   bash environments/run_all_validations.sh [--skip-slurm] [--quick]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Parse arguments
SKIP_SLURM=false
QUICK_MODE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-slurm)
            SKIP_SLURM=true
            shift
            ;;
        --quick)
            QUICK_MODE=true
            shift
            ;;
        --help)
            cat << EOF
Usage: $0 [OPTIONS]

Run all validation tests for Task 11 (Library import validation)

Options:
  --skip-slurm    Skip SLURM integration tests (for login node testing)
  --quick         Run quick tests only (skip comprehensive functionality tests)
  --help          Show this help message

Test Suites:
  11.1 - Python package imports (login node)
  11.2 - R package loading (login node)
  11.3 - SLURM integration (requires compute nodes)
  11.4 - Pipeline script compatibility

Examples:
  $0                    # Run all tests
  $0 --skip-slurm       # Run only login node tests
  $0 --quick            # Run quick validation only
EOF
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

print_header() {
    echo ""
    echo -e "${BLUE}======================================================================${NC}"
    echo -e "${BLUE}  $1${NC}"
    echo -e "${BLUE}======================================================================${NC}"
}

print_section() {
    echo ""
    echo -e "${BLUE}--- $1 ---${NC}"
}

print_success() {
    echo -e "${GREEN}✓${NC} $1"
}

print_failure() {
    echo -e "${RED}✗${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

# Test results tracking
declare -A TEST_RESULTS
TEST_ORDER=()

run_test() {
    local test_name="$1"
    local test_command="$2"
    
    TEST_ORDER+=("$test_name")
    
    print_section "$test_name"
    
    if eval "$test_command"; then
        TEST_RESULTS["$test_name"]="PASSED"
        print_success "$test_name completed successfully"
        return 0
    else
        TEST_RESULTS["$test_name"]="FAILED"
        print_failure "$test_name failed"
        return 1
    fi
}

# ============================================================================
# Task 11.1: Validate Python Package Imports
# ============================================================================
run_python_validation() {
    print_header "Task 11.1: Python Package Import Validation"
    
    if [ ! -f "${SCRIPT_DIR}/validate_python_imports.py" ]; then
        print_failure "Python validation script not found"
        return 1
    fi
    
    # Make executable
    chmod +x "${SCRIPT_DIR}/validate_python_imports.py"
    
    # Run on login node
    echo "Running Python validation on login node..."
    if ${SCRIPT_DIR}/run_with_env.sh python "${SCRIPT_DIR}/validate_python_imports.py"; then
        print_success "Python packages validated on login node"
        return 0
    else
        print_failure "Python validation failed on login node"
        return 1
    fi
}

# ============================================================================
# Task 11.2: Validate R Package Loading
# ============================================================================
run_r_validation() {
    print_header "Task 11.2: R Package Loading Validation"
    
    if [ ! -f "${SCRIPT_DIR}/validate_r_packages.R" ]; then
        print_failure "R validation script not found"
        return 1
    fi
    
    # Make executable
    chmod +x "${SCRIPT_DIR}/validate_r_packages.R"
    
    local all_passed=true
    
    # Test batch-effects environment
    echo "Running R validation on batch-effects environment (login node)..."
    if ${SCRIPT_DIR}/run_with_env.sh --r-env batch-effects Rscript "${SCRIPT_DIR}/validate_r_packages.R"; then
        print_success "R packages validated in batch-effects environment"
    else
        print_failure "R validation failed in batch-effects environment"
        all_passed=false
    fi
    
    # Test combatseq environment (if it exists)
    if [ -d "${SCRIPT_DIR}/r/combatseq" ]; then
        echo ""
        echo "Running R validation on combatseq environment (login node)..."
        if ${SCRIPT_DIR}/run_with_env.sh --r-env combatseq Rscript "${SCRIPT_DIR}/validate_r_packages.R"; then
            print_success "R packages validated in combatseq environment"
        else
            print_warning "R validation failed in combatseq environment (non-critical)"
        fi
    else
        print_warning "ComBat-seq environment not found, skipping"
    fi
    
    if [ "$all_passed" = true ]; then
        return 0
    else
        return 1
    fi
}

# ============================================================================
# Task 11.3: Validate SLURM Integration
# ============================================================================
run_slurm_validation() {
    print_header "Task 11.3: SLURM Integration Validation"
    
    if [ "$SKIP_SLURM" = true ]; then
        print_warning "Skipping SLURM tests (--skip-slurm flag)"
        return 0
    fi
    
    if ! command -v sbatch &>/dev/null; then
        print_warning "sbatch not available, skipping SLURM tests"
        return 0
    fi
    
    if [ ! -f "${SCRIPT_DIR}/validate_slurm_integration.sh" ]; then
        print_failure "SLURM validation script not found"
        return 1
    fi
    
    # Make executable
    chmod +x "${SCRIPT_DIR}/validate_slurm_integration.sh"
    
    # Run SLURM integration tests
    if bash "${SCRIPT_DIR}/validate_slurm_integration.sh"; then
        print_success "SLURM integration validated"
        return 0
    else
        print_failure "SLURM integration validation failed"
        return 1
    fi
}

# ============================================================================
# Task 11.4: Validate Pipeline Script Compatibility
# ============================================================================
run_pipeline_validation() {
    print_header "Task 11.4: Pipeline Script Compatibility Validation"
    
    if [ ! -f "${SCRIPT_DIR}/validate_pipeline_compatibility.sh" ]; then
        print_failure "Pipeline validation script not found"
        return 1
    fi
    
    # Make executable
    chmod +x "${SCRIPT_DIR}/validate_pipeline_compatibility.sh"
    
    # Run pipeline compatibility tests
    if bash "${SCRIPT_DIR}/validate_pipeline_compatibility.sh"; then
        print_success "Pipeline scripts validated"
        return 0
    else
        print_failure "Pipeline validation failed"
        return 1
    fi
}

# ============================================================================
# Main Execution
# ============================================================================
main() {
    local start_time=$(date +%s)
    
    print_header "Task 11: Library Import Validation - Master Test Suite"
    echo "  Started: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "  Mode: $([ "$QUICK_MODE" = true ] && echo "Quick" || echo "Comprehensive")"
    echo "  SLURM tests: $([ "$SKIP_SLURM" = true ] && echo "Disabled" || echo "Enabled")"
    echo ""
    
    # Check prerequisites
    if [ ! -f "${SCRIPT_DIR}/run_with_env.sh" ]; then
        echo "Error: run_with_env.sh not found at ${SCRIPT_DIR}/run_with_env.sh"
        exit 1
    fi
    
    if [ ! -f "${SCRIPT_DIR}/init_env.sh" ]; then
        echo "Error: init_env.sh not found at ${SCRIPT_DIR}/init_env.sh"
        exit 1
    fi
    
    # Source environment
    source "${SCRIPT_DIR}/init_env.sh"
    
    # Run test suites
    run_test "Task 11.1: Python Package Imports" "run_python_validation" || true
    run_test "Task 11.2: R Package Loading" "run_r_validation" || true
    
    if [ "$SKIP_SLURM" = false ]; then
        run_test "Task 11.3: SLURM Integration" "run_slurm_validation" || true
    fi
    
    run_test "Task 11.4: Pipeline Compatibility" "run_pipeline_validation" || true
    
    # Generate summary
    print_header "Final Summary"
    
    local total_tests=0
    local passed_tests=0
    local failed_tests=0
    
    echo ""
    echo "Test Results:"
    echo ""
    
    for test_name in "${TEST_ORDER[@]}"; do
        local result="${TEST_RESULTS[$test_name]}"
        total_tests=$((total_tests + 1))
        
        if [ "$result" = "PASSED" ]; then
            echo -e "  ${GREEN}✓${NC} $test_name"
            passed_tests=$((passed_tests + 1))
        else
            echo -e "  ${RED}✗${NC} $test_name"
            failed_tests=$((failed_tests + 1))
        fi
    done
    
    echo ""
    echo "Statistics:"
    echo "  Total tests:   $total_tests"
    echo "  Passed:        $passed_tests"
    echo "  Failed:        $failed_tests"
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    echo ""
    echo "  Duration:      ${duration}s"
    echo "  Completed:     $(date '+%Y-%m-%d %H:%M:%S')"
    echo ""
    
    if [ $failed_tests -eq 0 ]; then
        echo -e "${GREEN}✓ ALL VALIDATION TESTS PASSED${NC}"
        echo -e "${BLUE}======================================================================${NC}"
        echo ""
        echo "Task 11 (Library import validation) is complete!"
        echo ""
        echo "Next steps:"
        echo "  - Review test outputs for any warnings"
        echo "  - Proceed to Task 12 (Performance benchmarking)"
        echo "  - Update documentation with validation results"
        echo ""
        return 0
    else
        echo -e "${RED}✗ SOME VALIDATION TESTS FAILED${NC}"
        echo -e "${BLUE}======================================================================${NC}"
        echo ""
        echo "Please review the failed tests above and:"
        echo "  1. Check environment setup (init_env.sh)"
        echo "  2. Verify Python environment (uv sync)"
        echo "  3. Verify R environments (nix-shell)"
        echo "  4. Check SLURM configuration"
        echo ""
        return 1
    fi
}

# Run main function
main
exit $?
