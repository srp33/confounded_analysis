#!/bin/bash
# regenerate_minimal.sh
# Re-run Phase 1 with minimal package list and GitHub sources
# This script regenerates default.nix and .Rprofile with a minimal core set

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Regenerating batch-effects environment with minimal package list ==="
echo ""
echo "This will:"
echo "  1. Use generate_env_minimal.R (core packages only)"
echo "  2. Add GitHub sources for preprocessCore and polyester"
echo "  3. Regenerate default.nix and .Rprofile"
echo ""
echo "Press Enter to continue or Ctrl+C to cancel..."
read

# Nix configuration
NIX_ROOT="/grphome/grp_batch_effects/nix"
NIX_CHROOT_CMD="$NIX_ROOT/nix-user-chroot $NIX_ROOT"

echo "--- Entering Authoring Shell to generate default.nix ---"
echo ""

# Run Phase 1 with minimal package list
$NIX_CHROOT_CMD bash -c " \
    source ~/.nix-profile/etc/profile.d/nix.sh && \
    echo 'Running generate_env_minimal.R...' && \
    nix-shell \
      --option substituters 'https://cache.nixos.org https://rstats-on-nix.cachix.org' \
      --option trusted-public-keys 'cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= rstats-on-nix.cachix.org-1:9cJb6nqYZgKqgH5XJQN8FPkXqKlGqKqJqKqKqKqKqKo=' \
      -p R rPackages.rix --run 'Rscript generate_env_minimal.R' \
"

echo ""
echo "--- Generation complete ---"
echo ""
echo "Files created/updated:"
echo "  - default.nix (minimal package list)"
echo "  - .Rprofile (library isolation)"
echo ""
echo "Next steps:"
echo "  1. Review default.nix to verify package list"
echo "  2. Build the environment: ./build_with_cache.sh"
echo "  3. Install additional packages: nix-shell --run 'Rscript setup_packages.R'"
echo "  4. Test: nix-shell --run 'R --version'"
echo ""
echo "Note: This minimal environment includes ~30 core packages."
echo "      Additional packages will be installed via R's package manager."
