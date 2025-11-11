#!/bin/bash
# test_environment.sh
# Comprehensive test suite for the combatseq R environment

set -euo pipefail

# Configuration
NIX_ROOT="/grphome/grp_batch_effects/nix"
NIX_CHROOT_CMD="$NIX_ROOT/nix-user-chroot"
ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Testing combatseq R Environment ==="
echo "Environment directory: $ENV_DIR"
echo ""

# Track test results
TESTS_PASSED=0
TESTS_FAILED=0

# Function to run a test
run_test() {
    local test_name="$1"
    local test_cmd="$2"
    
    echo "--- Test: $test_name ---"
    if eval "$test_cmd"; then
        echo "✓ PASSED: $test_name"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        echo ""
        return 0
    else
        echo "✗ FAILED: $test_name"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        echo ""
        return 1
    fi
}

# Test 1: R version
run_test "R version check" \
    "$NIX_CHROOT_CMD $NIX_ROOT bash -c \"
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        cd $ENV_DIR && \
        nix-shell --run 'R --version'
    \""

# Test 2: Package loading test
run_test "Package loading test" \
    "$NIX_CHROOT_CMD $NIX_ROOT bash -c \"
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        cd $ENV_DIR && \
        nix-shell --run 'Rscript test_packages.R'
    \""

# Test 3: ComBat-seq functionality test
run_test "ComBat-seq functionality test" \
    "$NIX_CHROOT_CMD $NIX_ROOT bash -c \"
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        cd $ENV_DIR && \
        nix-shell --run 'Rscript test_combatseq.R'
    \""

# Test 4: Library path isolation
run_test "Library path isolation" \
    "$NIX_CHROOT_CMD $NIX_ROOT bash -c \"
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        cd $ENV_DIR && \
        nix-shell --run 'Rscript -e \".libPaths()\"'
    \" | grep -q '/nix/store'"

# Test 5: Bioconductor version check
run_test "Bioconductor version check" \
    "$NIX_CHROOT_CMD $NIX_ROOT bash -c \"
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        cd $ENV_DIR && \
        nix-shell --run 'Rscript -e \"library(BiocManager); cat(as.character(BiocManager::version()))\"'
    \""

# Test 6: Core ComBat-seq packages
run_test "Core ComBat-seq packages (sva, limma, DESeq2)" \
    "$NIX_CHROOT_CMD $NIX_ROOT bash -c \"
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        cd $ENV_DIR && \
        nix-shell --run 'Rscript -e \"library(sva); library(limma); library(DESeq2); cat(\\\"All packages loaded successfully\\\n\\\")\"'
    \""

# Test 7: Nix store size
echo "--- Test: Nix store size ---"
STORE_SIZE=$(du -sh "$NIX_ROOT/store" 2>/dev/null | cut -f1 || echo "unknown")
echo "Nix store size: $STORE_SIZE"
echo "✓ PASSED: Nix store size"
TESTS_PASSED=$((TESTS_PASSED + 1))
echo ""

# Summary
echo "=== Test Summary ==="
echo "Tests passed: $TESTS_PASSED"
echo "Tests failed: $TESTS_FAILED"
echo ""

if [ $TESTS_FAILED -eq 0 ]; then
    echo "✓ All tests passed!"
    echo "Environment is ready for ComBat-seq workflows"
    exit 0
else
    echo "✗ Some tests failed"
    echo "Please review the output above for details"
    exit 1
fi
