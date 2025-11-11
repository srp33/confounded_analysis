#!/bin/bash
# quick_test.sh
# Quick test to verify the environment is building/working

set -euo pipefail

NIX_ROOT="/grphome/grp_batch_effects/nix"
ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Quick Environment Test ==="
echo ""

# Test 1: Check if nix-user-chroot works
echo "Test 1: Checking nix-user-chroot..."
if [ -f "$NIX_ROOT/nix-user-chroot" ]; then
    echo "✓ nix-user-chroot found"
else
    echo "✗ nix-user-chroot not found"
    exit 1
fi

# Test 2: Check if default.nix exists
echo "Test 2: Checking default.nix..."
if [ -f "$ENV_DIR/default.nix" ]; then
    echo "✓ default.nix found"
else
    echo "✗ default.nix not found"
    exit 1
fi

# Test 3: Check if .Rprofile exists
echo "Test 3: Checking .Rprofile..."
if [ -f "$ENV_DIR/.Rprofile" ]; then
    echo "✓ .Rprofile found"
else
    echo "✗ .Rprofile not found"
fi

# Test 4: Try to parse default.nix (syntax check)
echo "Test 4: Checking default.nix syntax..."
cd "$ENV_DIR"
if "$NIX_ROOT/nix-user-chroot" "$NIX_ROOT" bash -c "
    source ~/.nix-profile/etc/profile.d/nix.sh && \
    cd $ENV_DIR && \
    nix-instantiate --parse default.nix > /dev/null 2>&1
"; then
    echo "✓ default.nix syntax is valid"
else
    echo "✗ default.nix has syntax errors"
    exit 1
fi

# Test 5: Check if environment can be instantiated (doesn't build, just checks)
echo "Test 5: Checking if environment can be instantiated..."
cd "$ENV_DIR"
if "$NIX_ROOT/nix-user-chroot" "$NIX_ROOT" bash -c "
    source ~/.nix-profile/etc/profile.d/nix.sh && \
    cd $ENV_DIR && \
    nix-instantiate default.nix > /dev/null 2>&1
"; then
    echo "✓ Environment can be instantiated"
else
    echo "✗ Environment instantiation failed"
    exit 1
fi

echo ""
echo "=== All quick tests passed ==="
echo ""
echo "The environment is ready to build."
echo "To build and test the full environment, run:"
echo "  ./build_environment.sh"
echo ""
echo "This will take 20-40 minutes for the first build."
