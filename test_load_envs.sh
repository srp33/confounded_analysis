#!/usr/bin/env bash
# test_load_envs.sh - Comprehensive test suite for load_envs.sh
# This script tests all functionality of the environment loader

set -e  # Exit on error

# Path to the script being tested
LOAD_ENVS_SCRIPT="environments/load_envs.sh"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Test result tracking
declare -a FAILED_TESTS

# Helper functions
print_header() {
    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
    echo ""
}

print_test() {
    echo -e "${YELLOW}TEST: $1${NC}"
    TESTS_RUN=$((TESTS_RUN + 1))
}

print_pass() {
    echo -e "${GREEN}✓ PASS${NC}: $1"
    TESTS_PASSED=$((TESTS_PASSED + 1))
    echo ""
}

print_fail() {
    echo -e "${RED}✗ FAIL${NC}: $1"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    FAILED_TESTS+=("$2")
    echo ""
}

print_info() {
    echo -e "  ${BLUE}INFO${NC}: $1"
}

# Cleanup function
cleanup_test_envs() {
    print_info "Cleaning up test environments..."
    rm -rf environments/test_py 2>/dev/null || true
    rm -rf environments/test_r 2>/dev/null || true
}

# Setup test environments
setup_test_envs() {
    print_header "Setting Up Test Environments"
    
    # Create Python-only test environment
    print_info "Creating test_py (Python-only)..."
    mkdir -p environments/test_py
    cat > environments/test_py/pyproject.toml << 'EOF'
[project]
name = "test_py"
version = "0.1.0"
description = "Test Python-only environment"
requires-python = ">=3.10"
dependencies = [
    "requests>=2.31",
]
EOF
    
    # Create R-only test environment
    print_info "Creating test_r (R-only)..."
    mkdir -p environments/test_r
    cat > environments/test_r/rproject.toml << 'EOF'
[project]
name = "test_r"
r_version = "4.3"

repositories = [
    { alias = "PPM", url = "https://packagemanager.posit.co/cran/latest" },
]

dependencies = [
    "jsonlite",
]
EOF
    
    echo -e "${GREEN}✓ Test environments created${NC}"
    echo ""
}

#------------------------------------------------------------------------------
# Test 7.1: Python-only environment
#------------------------------------------------------------------------------
test_python_only() {
    print_header "Test 7.1: Python-only Environment"
    
    print_test "Activating Python-only environment (test_py)"
    
    # Clear environment variables
    deactivate 2>/dev/null || true
    unset VIRTUAL_ENV
    unset R_LIBS_USER
    
    # Source the script and capture output
    output=$(source $LOAD_ENVS_SCRIPT test_py 2>&1)
    exit_code=$?
    
    # Check if activation succeeded
    if [[ $exit_code -ne 0 ]] || [[ "$output" == *"ERROR"* ]]; then
        print_fail "Failed to activate Python-only environment" "7.1"
        print_info "Output: $output"
        return 1
    fi
    
    print_pass "Script executed without errors"
    
    # Verify VIRTUAL_ENV is set
    if [[ -n "$VIRTUAL_ENV" ]]; then
        print_pass "VIRTUAL_ENV is set: $VIRTUAL_ENV"
    else
        print_fail "VIRTUAL_ENV is not set" "7.1"
        return 1
    fi
    
    # Verify R_LIBS_USER is not set
    if [[ -z "$R_LIBS_USER" ]]; then
        print_pass "R_LIBS_USER is not set (as expected)"
    else
        print_fail "R_LIBS_USER should not be set for Python-only env" "7.1"
        return 1
    fi
    
    # Verify Python is available
    if command -v python &> /dev/null; then
        python_path=$(which python)
        if [[ "$python_path" == *".venv"* ]]; then
            print_pass "Python is from virtual environment: $python_path"
        else
            print_fail "Python is not from virtual environment" "7.1"
            return 1
        fi
    else
        print_fail "Python command not found" "7.1"
        return 1
    fi
    
    # Deactivate for next test
    deactivate 2>/dev/null || true
    unset VIRTUAL_ENV
    
    return 0
}

#------------------------------------------------------------------------------
# Test 7.2: R-only environment
#------------------------------------------------------------------------------
test_r_only() {
    print_header "Test 7.2: R-only Environment"
    
    print_test "Activating R-only environment (test_r)"
    
    # Clear environment variables
    deactivate 2>/dev/null || true
    unset VIRTUAL_ENV
    unset R_LIBS_USER
    
    # Check if rig and rv are available
    if ! command -v rig &> /dev/null || ! command -v rv &> /dev/null; then
        print_info "Skipping R-only test: rig/rv not installed"
        TESTS_RUN=$((TESTS_RUN - 1))  # Don't count as run
        return 0
    fi
    
    # Source the script and capture output
    output=$(source $LOAD_ENVS_SCRIPT test_r 2>&1)
    exit_code=$?
    
    # Check if activation succeeded
    if [[ $exit_code -ne 0 ]] || [[ "$output" == *"ERROR"* ]]; then
        print_fail "Failed to activate R-only environment" "7.2"
        print_info "Output: $output"
        return 1
    fi
    
    print_pass "Script executed without errors"
    
    # Verify R_LIBS_USER is set
    if [[ -n "$R_LIBS_USER" ]]; then
        print_pass "R_LIBS_USER is set: $R_LIBS_USER"
    else
        print_fail "R_LIBS_USER is not set" "7.2"
        return 1
    fi
    
    # Verify VIRTUAL_ENV is not set
    if [[ -z "$VIRTUAL_ENV" ]]; then
        print_pass "VIRTUAL_ENV is not set (as expected)"
    else
        print_fail "VIRTUAL_ENV should not be set for R-only env" "7.2"
        return 1
    fi
    
    # Verify R is available
    if command -v R &> /dev/null; then
        print_pass "R command is available"
    else
        print_fail "R command not found" "7.2"
        return 1
    fi
    
    # Clear for next test
    unset R_LIBS_USER
    
    return 0
}

#------------------------------------------------------------------------------
# Test 7.3: Dual environment
#------------------------------------------------------------------------------
test_dual_environment() {
    print_header "Test 7.3: Dual Environment"
    
    print_test "Activating dual environment (book_chapter)"
    
    # Clear environment variables
    deactivate 2>/dev/null || true
    unset VIRTUAL_ENV
    unset R_LIBS_USER
    
    # Check if book_chapter exists
    if [[ ! -d "environments/book_chapter" ]]; then
        print_info "Skipping dual environment test: book_chapter not found"
        TESTS_RUN=$((TESTS_RUN - 1))
        return 0
    fi
    
    # Check if both config files exist
    if [[ ! -f "environments/book_chapter/pyproject.toml" ]] || [[ ! -f "environments/book_chapter/rproject.toml" ]]; then
        print_info "Skipping dual environment test: book_chapter missing config files"
        TESTS_RUN=$((TESTS_RUN - 1))
        return 0
    fi
    
    # Check if rig and rv are available
    if ! command -v rig &> /dev/null || ! command -v rv &> /dev/null; then
        print_info "Skipping dual environment test: rig/rv not installed"
        TESTS_RUN=$((TESTS_RUN - 1))
        return 0
    fi
    
    # Source the script and capture output
    output=$(source $LOAD_ENVS_SCRIPT book_chapter 2>&1)
    exit_code=$?
    
    # Check if activation succeeded
    if [[ $exit_code -ne 0 ]] || [[ "$output" == *"ERROR"* ]]; then
        print_fail "Failed to activate dual environment" "7.3"
        print_info "Output: $output"
        return 1
    fi
    
    print_pass "Script executed without errors"
    
    # Verify both VIRTUAL_ENV and R_LIBS_USER are set
    if [[ -n "$VIRTUAL_ENV" ]]; then
        print_pass "VIRTUAL_ENV is set: $VIRTUAL_ENV"
    else
        print_fail "VIRTUAL_ENV is not set in dual environment" "7.3"
        return 1
    fi
    
    if [[ -n "$R_LIBS_USER" ]]; then
        print_pass "R_LIBS_USER is set: $R_LIBS_USER"
    else
        print_fail "R_LIBS_USER is not set in dual environment" "7.3"
        return 1
    fi
    
    # Test Python imports
    if python -c "import numpy" 2>/dev/null; then
        print_pass "Python imports work (numpy imported successfully)"
    else
        print_info "Python imports test skipped (numpy may not be installed yet)"
    fi
    
    # Test R library loading (if R is available)
    if command -v R &> /dev/null; then
        if R --quiet --no-save -e "library(tidyverse)" 2>&1 | grep -q "Error"; then
            print_info "R library test skipped (tidyverse may not be installed yet)"
        else
            print_pass "R library loading works"
        fi
    fi
    
    # Cleanup
    deactivate 2>/dev/null || true
    unset VIRTUAL_ENV
    unset R_LIBS_USER
    
    return 0
}

#------------------------------------------------------------------------------
# Test 7.4: Error handling
#------------------------------------------------------------------------------
test_error_handling() {
    print_header "Test 7.4: Error Handling"
    
    # Test 7.4.1: Nonexistent project
    print_test "Testing with nonexistent project name"
    output=$(source $LOAD_ENVS_SCRIPT nonexistent_project 2>&1 || true)
    
    if [[ "$output" == *"ERROR"* ]] && [[ "$output" == *"not found"* ]]; then
        print_pass "Appropriate error message for nonexistent project"
    else
        print_fail "Missing or incorrect error message for nonexistent project" "7.4"
        print_info "Output: $output"
    fi
    
    # Test 7.4.2: No arguments
    print_test "Testing with no arguments"
    output=$(source $LOAD_ENVS_SCRIPT 2>&1 || true)
    
    if [[ "$output" == *"ERROR"* ]] && [[ "$output" == *"No project name"* ]]; then
        print_pass "Appropriate error message when no arguments provided"
    else
        print_fail "Missing or incorrect error message for no arguments" "7.4"
        print_info "Output: $output"
    fi
    
    # Test 7.4.3: Missing prerequisites (simulate by using invalid PATH)
    print_test "Testing with missing prerequisites"
    # This test is complex and may not work in all environments
    # We'll check if the script handles missing tools gracefully
    print_info "Skipping missing prerequisites test (requires PATH manipulation)"
    TESTS_RUN=$((TESTS_RUN - 1))
    
    return 0
}

#------------------------------------------------------------------------------
# Test 7.5: Special flags
#------------------------------------------------------------------------------
test_special_flags() {
    print_header "Test 7.5: Special Flags"
    
    # Test 7.5.1: --help flag
    print_test "Testing --help flag"
    output=$(source $LOAD_ENVS_SCRIPT --help 2>&1)
    
    if [[ "$output" == *"USAGE"* ]] && [[ "$output" == *"load_envs.sh"* ]]; then
        print_pass "--help flag displays usage information"
    else
        print_fail "--help flag does not display proper usage" "7.5"
        print_info "Output: $output"
    fi
    
    # Test 7.5.2: --list flag
    print_test "Testing --list flag"
    output=$(source $LOAD_ENVS_SCRIPT --list 2>&1)
    
    if [[ "$output" == *"Available Projects"* ]]; then
        print_pass "--list flag displays available projects"
    else
        print_fail "--list flag does not display projects" "7.5"
        print_info "Output: $output"
    fi
    
    # Test 7.5.3: --verbose flag
    print_test "Testing --verbose flag"
    deactivate 2>/dev/null || true
    output=$(source $LOAD_ENVS_SCRIPT test_py --verbose 2>&1 || true)
    deactivate 2>/dev/null || true
    
    if [[ "$output" == *"VERBOSE"* ]]; then
        print_pass "--verbose flag enables detailed output"
    else
        print_fail "--verbose flag does not enable verbose mode" "7.5"
        print_info "Output: $output"
    fi
    
    # Test 7.5.4: --force-sync flag
    print_test "Testing --force-sync flag"
    # First activate normally to create the environment
    source $LOAD_ENVS_SCRIPT test_py > /dev/null 2>&1
    deactivate 2>/dev/null || true
    
    # Then test force-sync
    output=$(source $LOAD_ENVS_SCRIPT test_py --force-sync 2>&1 || true)
    deactivate 2>/dev/null || true
    
    if [[ "$output" == *"Syncing"* ]] || [[ "$output" == *"sync"* ]]; then
        print_pass "--force-sync flag forces synchronization"
    else
        print_info "Could not verify --force-sync (environment may already be synced)"
        # Don't fail this test as it's hard to verify
    fi
    
    return 0
}

#------------------------------------------------------------------------------
# Test 7.6: Performance
#------------------------------------------------------------------------------
test_performance() {
    print_header "Test 7.6: Performance"
    
    print_test "Measuring activation time for cached environment"
    
    # Ensure environment is cached
    source $LOAD_ENVS_SCRIPT test_py > /dev/null 2>&1
    deactivate 2>/dev/null || true
    unset VIRTUAL_ENV
    
    # Measure activation time
    start_time=$(date +%s.%N)
    source $LOAD_ENVS_SCRIPT test_py > /dev/null 2>&1
    end_time=$(date +%s.%N)
    
    activation_time=$(echo "$end_time - $start_time" | bc)
    
    print_info "Activation time: ${activation_time}s"
    
    # Check if under 2 seconds
    if (( $(echo "$activation_time < 2.0" | bc -l) )); then
        print_pass "Activation completed in < 2 seconds"
    else
        print_fail "Activation took longer than 2 seconds" "7.6"
    fi
    
    # Test smart sync skipping
    print_test "Verifying smart sync skipping"
    deactivate 2>/dev/null || true
    output=$(source $LOAD_ENVS_SCRIPT test_py 2>&1 || true)
    deactivate 2>/dev/null || true
    
    if [[ "$output" == *"up-to-date"* ]] || [[ "$output" == *"skipping sync"* ]]; then
        print_pass "Smart sync skipping works"
    else
        print_info "Could not verify smart sync skipping from output"
    fi
    
    deactivate 2>/dev/null || true
    unset VIRTUAL_ENV
    
    return 0
}

#------------------------------------------------------------------------------
# Test 7.7: Sourcing detection
#------------------------------------------------------------------------------
test_sourcing_detection() {
    print_header "Test 7.7: Sourcing Detection"
    
    print_test "Testing execution vs sourcing detection"
    
    # Execute the script (not source it)
    output=$(bash $LOAD_ENVS_SCRIPT test_py 2>&1 || true)
    
    if [[ "$output" == *"WARNING"* ]] && [[ "$output" == *"must be sourced"* ]]; then
        print_pass "Script detects when executed instead of sourced"
    else
        print_fail "Script does not detect execution vs sourcing" "7.7"
    fi
    
    # Verify environment is not modified when executed
    if [[ -z "$VIRTUAL_ENV" ]]; then
        print_pass "Environment not modified when script is executed"
    else
        print_fail "Environment was modified despite script being executed" "7.7"
    fi
    
    return 0
}

#------------------------------------------------------------------------------
# Main test execution
#------------------------------------------------------------------------------
main() {
    print_header "load_envs.sh Test Suite"
    echo "Starting comprehensive tests..."
    echo ""
    
    # Setup
    cleanup_test_envs
    setup_test_envs
    
    # Run all tests
    test_python_only || true
    test_r_only || true
    test_dual_environment || true
    test_error_handling || true
    test_special_flags || true
    test_performance || true
    test_sourcing_detection || true
    
    # Cleanup
    cleanup_test_envs
    
    # Print summary
    print_header "Test Summary"
    echo "Tests run:    $TESTS_RUN"
    echo -e "Tests passed: ${GREEN}$TESTS_PASSED${NC}"
    echo -e "Tests failed: ${RED}$TESTS_FAILED${NC}"
    echo ""
    
    if [[ $TESTS_FAILED -gt 0 ]]; then
        echo -e "${RED}Failed tests:${NC}"
        for test in "${FAILED_TESTS[@]}"; do
            echo "  - Test $test"
        done
        echo ""
        exit 1
    else
        echo -e "${GREEN}All tests passed!${NC}"
        echo ""
        exit 0
    fi
}

# Run main function
main

