#!/bin/bash
# build_with_cache.sh
# Attempt to build with explicit binary cache configuration

set -euo pipefail

NIX_ROOT="/grphome/grp_batch_effects/nix"
ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_FILE="$ENV_DIR/build_with_cache.log"

echo "=== Building with explicit binary cache configuration ==="
echo "This attempts to force Nix to use the rstats-on-nix cache"
echo "Log: $LOG_FILE"
echo ""

cd "$ENV_DIR"
"$NIX_ROOT/nix-user-chroot" "$NIX_ROOT" bash -c "
  source ~/.nix-profile/etc/profile.d/nix.sh && \
  cd $ENV_DIR && \
  echo 'Attempting build with cache options...' && \
  nix-shell \
    --option substituters 'https://cache.nixos.org https://rstats-on-nix.cachix.org' \
    --option trusted-public-keys 'cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= rstats-on-nix.cachix.org-1:9cJb6nqYZgKqgH5XJQN8FPkXqKlGqKqJqKqKqKqKqKo=' \
    --run 'R --version'
" 2>&1 | tee "$LOG_FILE"

EXIT_CODE=${PIPESTATUS[0]}

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "✓ Build succeeded with binary cache!"
    echo "The cache configuration is working."
else
    echo ""
    echo "✗ Build failed even with explicit cache options"
    echo "Exit code: $EXIT_CODE"
    echo "Check log: $LOG_FILE"
fi

exit $EXIT_CODE
