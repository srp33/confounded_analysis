#!/bin/bash
# nix_activate.sh
# Phase 2 (Activation): Use the generated default.nix to activate the R environment
# This script provides an interactive shell with the batch-effects R environment

set -euo pipefail

# Configuration
NIX_ROOT="/grphome/grp_batch_effects/nix"
NIX_CHROOT_CMD="$NIX_ROOT/nix-user-chroot"
ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Activating batch-effects R environment (Phase 2) ==="
echo "Environment directory: $ENV_DIR"
echo ""

# Check if nix-user-chroot exists
if [ ! -f "$NIX_CHROOT_CMD" ]; then
    echo "ERROR: nix-user-chroot not found at $NIX_CHROOT_CMD"
    echo "Please ensure Nix is installed in the group directory."
    exit 1
fi

# Check if default.nix exists
if [ ! -f "$ENV_DIR/default.nix" ]; then
    echo "ERROR: default.nix not found in $ENV_DIR"
    echo "Please run Phase 1 first: ./run_generator.sh"
    exit 1
fi

# Check if .Rprofile exists
if [ ! -f "$ENV_DIR/.Rprofile" ]; then
    echo "WARNING: .Rprofile not found. Library isolation may not work correctly."
fi

echo "Starting nix-shell..."
echo "This will:"
echo "  1. Enter the nix-user-chroot namespace"
echo "  2. Source the Nix profile"
echo "  3. Activate the R environment from default.nix"
echo "  4. Load .Rprofile for library path isolation"
echo ""
echo "First-time build may take 20-40 minutes to download packages."
echo "Subsequent activations will be fast (~500ms)."
echo ""

# Enter nix-user-chroot and start nix-shell
# The --rcfile creates an interactive bash shell with nix-shell activated
cd "$ENV_DIR"
exec "$NIX_CHROOT_CMD" "$NIX_ROOT" bash --rcfile <(cat <<'EOF'
# Source the Nix profile to make nix-shell available
if [ -f ~/.nix-profile/etc/profile.d/nix.sh ]; then
    source ~/.nix-profile/etc/profile.d/nix.sh
else
    echo "ERROR: Nix profile not found at ~/.nix-profile/etc/profile.d/nix.sh"
    exit 1
fi

# Activate nix-shell (this loads default.nix)
echo "Entering nix-shell..."
nix-shell

# If nix-shell exits, return to normal bash
echo ""
echo "Exited nix-shell. Type 'exit' to leave the nix-user-chroot namespace."
EOF
)
