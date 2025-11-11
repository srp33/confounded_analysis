#!/bin/bash
# build_environment.sh
# Build the combatseq R environment for the first time
# This script initiates the Nix build process and logs progress

set -euo pipefail

# Configuration
NIX_ROOT="/grphome/grp_batch_effects/nix"
NIX_CHROOT_CMD="$NIX_ROOT/nix-user-chroot"
ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_FILE="$ENV_DIR/build.log"

echo "=== Building combatseq R Environment ==="
echo "Environment directory: $ENV_DIR"
echo "Log file: $LOG_FILE"
echo ""

# Check prerequisites
if [ ! -f "$NIX_CHROOT_CMD" ]; then
    echo "ERROR: nix-user-chroot not found at $NIX_CHROOT_CMD"
    exit 1
fi

if [ ! -f "$ENV_DIR/default.nix" ]; then
    echo "ERROR: default.nix not found in $ENV_DIR"
    exit 1
fi

echo "Starting build process..."
echo "This will take 20-40 minutes for the first build."
echo "Progress will be logged to: $LOG_FILE"
echo ""
echo "You can monitor progress in another terminal with:"
echo "  tail -f $LOG_FILE"
echo ""

# Record start time
START_TIME=$(date +%s)

# Run the build
cd "$ENV_DIR"
"$NIX_CHROOT_CMD" "$NIX_ROOT" bash -c "
  source ~/.nix-profile/etc/profile.d/nix.sh && \
  cd $ENV_DIR && \
  nix-shell --run 'R --version'
" 2>&1 | tee "$LOG_FILE"

# Record end time
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

# Check if build succeeded
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo ""
    echo "=== Build Successful ==="
    echo "Duration: ${MINUTES}m ${SECONDS}s"
    echo "Log saved to: $LOG_FILE"
    echo ""
    echo "Next steps:"
    echo "  1. Run tests: ./test_environment.sh"
    echo "  2. Activate environment: ./nix_activate.sh"
    exit 0
else
    echo ""
    echo "=== Build Failed ==="
    echo "Duration: ${MINUTES}m ${SECONDS}s"
    echo "Check log file for details: $LOG_FILE"
    exit 1
fi
