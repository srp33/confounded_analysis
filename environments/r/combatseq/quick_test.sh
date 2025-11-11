#!/bin/bash
# quick_test.sh
# Quick validation of combatseq environment without building
# This script checks prerequisites and validates Nix syntax

set -euo pipefail

# Configuration
NIX_ROOT="/grphome/grp_batch_effects/nix"
NIX_CHROOT_CMD="$NIX_ROOT/nix-user-chroot"
ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Quick Test: combatseq Environment ==="
echo "Environment directory: $ENV_DIR"
echo ""

# Track test results
TESTS_PASSED=0
TESTS_FAILED=0

# Function to run a test
run_test() {
    local test_name="$1"
    local test_cmd="$2"
    
    echo -n "Testing $test_name... "
    if eval "$test_cmd" > /dev/null 2>&1; then
        echo "✓ PASSED"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo "✗ FAILED"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# Test 1: nix-user-chroot exists
run_test "nix-user-chroot availability" \
    "[ -f $NIX_CHROOT_CMD ]"

# Test 2: default.nix exists
run_test "default.nix existence" \
    "[ -f $ENV_DIR/default.nix ]"

# Test 3: .Rprofile exists
run_test ".Rprofile existence" \
    "[ -f $ENV_DIR/.Rprofile ]"

# Test 4: Nix syntax validation
run_test "default.nix syntax" \
    "$NIX_CHROOT_CMD $NIX_ROOT bash -c \"
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        nix-instantiate --parse $ENV_DIR/default.nix
    \""

# Test 5: Environment instantiation
echo -n "Testing environment instantiation... "
DERIVATIONS=$("$NIX_CHROOT_CMD" "$NIX_ROOT" bash -c "
    source ~/.nix-profile/etc/profile.d/nix.sh && \
    cd $ENV_DIR && \
    nix-instantiate --show-trace default.nix 2>&1
" | wc -l)

if [ $? -eq 0 ]; then
    echo "✓ PASSED ($DERIVATIONS derivations)"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    echo "✗ FAILED"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

# Summary
echo ""
echo "=== Quick Test Summary ==="
echo "Tests passed: $TESTS_PASSED"
echo "Tests failed: $TESTS_FAILED"
echo ""

if [ $TESTS_FAILED -eq 0 ]; then
    echo "✓ All quick tests passed!"
    echo "Environment is ready to build"
    echo ""
    echo "Next step: ./build_environment.sh"
    exit 0
else
    echo "✗ Some tests failed"
    echo "Please fix the issues before building"
    exit 1
fi
