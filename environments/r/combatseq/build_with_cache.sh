#!/bin/bash
# build_with_cache.sh
# Build the combatseq R environment with explicit binary cache configuration
# This ensures Nix downloads pre-built packages instead of building from source

set -euo pipefail

# Configuration
NIX_ROOT="/grphome/grp_batch_effects/nix"
NIX_CHROOT_CMD="$NIX_ROOT/nix-user-chroot"
ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_FILE="$ENV_DIR/build_with_cache.log"

echo "=== Building combatseq R Environment with Binary Cache ==="
echo "Environment directory: $ENV_DIR"
echo "Log file: $LOG_FILE"
echo ""

# Check prerequisites
if [ ! -f "$NIX_CHROOT_CMD" ]; then
    echo "ERROR: nix-user-chroot not found at $NIX_CHROOT_CMD"
    echo "Please ensure Nix is installed in the group directory."
    exit 1
fi

if [ ! -f "$ENV_DIR/default.nix" ]; then
    echo "ERROR: default.nix not found in $ENV_DIR"
    echo "Please run Phase 1 first: ./run_generator.sh"
    exit 1
fi

echo "Starting build with binary cache configuration..."
echo "Cache: https://rstats-on-nix.cachix.org"
echo ""
echo "This will download pre-built packages (~10-20 minutes)"
echo "instead of building from source (~20-40 minutes)."
echo ""
echo "Progress will be logged to: $LOG_FILE"
echo ""
echo "You can monitor progress in another terminal with:"
echo "  tail -f $LOG_FILE"
echo ""

# Record start time
START_TIME=$(date +%s)

# Run the build with explicit cache options
cd "$ENV_DIR"
"$NIX_CHROOT_CMD" "$NIX_ROOT" bash -c "
  source ~/.nix-profile/etc/profile.d/nix.sh && \
  cd $ENV_DIR && \
  echo 'Building with cache options...' && \
  nix-shell \
    --option substituters 'https://cache.nixos.org https://rstats-on-nix.cachix.org' \
    --option trusted-public-keys 'cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= rstats-on-nix.cachix.org-1:vdiiVgocg6WeJrODIqdprZRUrhi1JzhBnXv7aWI6+F0=' \
    --run 'R --version'
" 2>&1 | tee "$LOG_FILE"

EXIT_CODE=${PIPESTATUS[0]}

# Record end time
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

# Check if build succeeded
if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "=== Build Successful ==="
    echo "Duration: ${MINUTES}m ${SECONDS}s"
    echo "Log saved to: $LOG_FILE"
    echo ""
    echo "The binary cache configuration is working!"
    echo ""
    echo "Next steps:"
    echo "  1. Test package loading: ./test_packages.R"
    echo "  2. Test ComBat-seq scripts: ./test_combatseq.R"
    echo "  3. Activate environment: ./nix_activate.sh"
    exit 0
else
    echo ""
    echo "=== Build Failed ==="
    echo "Duration: ${MINUTES}m ${SECONDS}s"
    echo "Exit code: $EXIT_CODE"
    echo "Check log file for details: $LOG_FILE"
    exit 1
fi
