#!/bin/bash
# build_nix_env.sh
# Pre-build Nix environment to avoid 5-minute activation overhead
#
# This script runs nix-build once to create a ./result symlink.
# Subsequent nix-shell invocations on ./result are nearly instant.
#
# Usage:
#   cd environments/r/combatseq  # or batch-effects
#   ../../build_nix_env.sh

set -euo pipefail

# Nix configuration
NIX_ROOT="/grphome/grp_batch_effects/nix"
NIX_CHROOT_CMD="$NIX_ROOT/nix-user-chroot $NIX_ROOT"

# Binary cache configuration - get from user's nix.conf
if [ -f ~/.config/nix/nix.conf ]; then
    echo "Using binary cache configuration from ~/.config/nix/nix.conf"
    # Extract the actual key from config
    CACHE_KEY=$(grep "rstats-on-nix.cachix.org-1:" ~/.config/nix/nix.conf | sed 's/.*rstats-on-nix.cachix.org-1:\([^ ]*\).*/\1/' || echo "")
    if [ -z "$CACHE_KEY" ]; then
        echo "⚠️  Warning: rstats-on-nix cache key not found in nix.conf"
        echo "Run: bash environments/r/fix_nix_cache.sh"
        exit 1
    fi
    NIX_CACHE_OPTS="--option substituters 'https://cache.nixos.org https://rstats-on-nix.cachix.org' --option trusted-public-keys 'cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= rstats-on-nix.cachix.org-1:${CACHE_KEY}'"
else
    echo "❌ Error: ~/.config/nix/nix.conf not found"
    echo "Run: bash environments/r/fix_nix_cache.sh"
    exit 1
fi

# Check we're in an R environment directory
if [ ! -f "default.nix" ]; then
    echo "❌ Error: default.nix not found in current directory" >&2
    echo "Run this script from an R environment directory (e.g., environments/r/combatseq)" >&2
    exit 1
fi

# Check nix-user-chroot exists
if [ ! -f "$NIX_ROOT/nix-user-chroot" ]; then
    echo "❌ Error: nix-user-chroot not found at $NIX_ROOT/nix-user-chroot" >&2
    echo "Nix installation may be incomplete" >&2
    exit 1
fi

ENV_DIR=$(pwd)
ENV_NAME=$(basename "$ENV_DIR")

echo "============================================"
echo "Building Nix environment: $ENV_NAME"
echo "Directory: $ENV_DIR"
echo "============================================"
echo ""

# Check if result already exists
if [ -L "result" ]; then
    echo "Found existing ./result symlink"
    echo "Target: $(readlink result)"
    echo ""
    read -p "Rebuild? This will take ~5 minutes. [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Skipping rebuild. Using existing result."
        exit 0
    fi
    echo ""
fi

echo "Running nix-build (this will take ~5 minutes on first run)..."
echo "Command: nix-build $NIX_CACHE_OPTS -A shell"
echo ""

START_TIME=$(date +%s)

# Run nix-build to create the result symlink
echo "Entering nix-user-chroot environment..."
if ! $NIX_CHROOT_CMD bash -c "
    set -euo pipefail
    source ~/.nix-profile/etc/profile.d/nix.sh || { echo '❌ Failed to source nix profile'; exit 1; }
    cd '$ENV_DIR' || { echo '❌ Failed to cd to $ENV_DIR'; exit 1; }
    echo 'Running nix-build...'
    nix-build $NIX_CACHE_OPTS -A shell 2>&1 | tee build.log
"; then
    echo ""
    echo "❌ Build failed!"
    echo ""
    echo "Check build.log for details"
    echo ""
    echo "Common issues:"
    echo "  1. Binary cache trust not configured - run: bash environments/r/fix_nix_cache.sh"
    echo "  2. Network issues - check connectivity to rstats-on-nix.cachix.org"
    echo "  3. Disk space - check available space in /grphome/grp_batch_effects/nix/"
    exit 1
fi

END_TIME=$(date +%s)
BUILD_TIME=$((END_TIME - START_TIME))

echo ""
echo "============================================"
echo "Build complete in ${BUILD_TIME}s"
echo "============================================"
echo ""
echo "Result symlink created:"
ls -lh result
echo ""
echo "Target store path:"
readlink result
echo ""
echo "Next steps:"
echo "  1. Test activation speed:"
echo "     time nix-shell result --run 'R --version'"
echo ""
echo "  2. Use with run_with_env.sh (automatic)"
echo "     The wrapper will detect ./result and use it"
echo ""
echo "Rebuild only when default.nix changes."
