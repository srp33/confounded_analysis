#!/bin/bash
# Configure Nix binary cache using cachix
# This is the recommended approach per rix documentation

set -euo pipefail

echo "=========================================="
echo "Nix Binary Cache Setup with cachix"
echo "=========================================="
echo ""

NIX_ROOT="/grphome/grp_batch_effects/nix"
NIX_CHROOT_CMD="$NIX_ROOT/nix-user-chroot $NIX_ROOT"

# Check if cachix is installed
echo "Checking for cachix..."
if ! $NIX_CHROOT_CMD bash -c "source ~/.nix-profile/etc/profile.d/nix.sh && command -v cachix" &>/dev/null; then
    echo "cachix not found. Installing..."
    echo ""
    
    $NIX_CHROOT_CMD bash -c "
        source ~/.nix-profile/etc/profile.d/nix.sh
        nix-env -iA cachix -f https://cachix.org/api/v1/install
    "
    
    if [ $? -eq 0 ]; then
        echo "✅ cachix installed successfully"
    else
        echo "❌ Failed to install cachix"
        exit 1
    fi
else
    echo "✅ cachix already installed"
fi

echo ""
echo "Configuring rstats-on-nix cache..."
echo ""

# Use cachix to configure the cache
# This is the authoritative method per rix documentation
$NIX_CHROOT_CMD bash -c "
    source ~/.nix-profile/etc/profile.d/nix.sh
    cachix use rstats-on-nix
"

if [ $? -eq 0 ]; then
    echo ""
    echo "✅ rstats-on-nix cache configured successfully!"
    echo ""
    echo "=========================================="
    echo "Verification"
    echo "=========================================="
    echo ""
    
    if [ -f ~/.config/nix/nix.conf ]; then
        echo "Configuration file: ~/.config/nix/nix.conf"
        echo ""
        echo "Substituters:"
        grep "substituters" ~/.config/nix/nix.conf || echo "  (using defaults)"
        echo ""
        echo "Trusted public keys:"
        grep "trusted-public-keys" ~/.config/nix/nix.conf || echo "  (using defaults)"
    else
        echo "Note: cachix may have configured system-level settings"
        echo "This is normal and expected"
    fi
    
    echo ""
    echo "=========================================="
    echo "Next Steps"
    echo "=========================================="
    echo ""
    echo "1. Build the R environment:"
    echo "   cd environments/r/batch-effects"
    echo "   bash ../build_nix_env.sh"
    echo ""
    echo "2. Expected behavior:"
    echo "   ✅ Packages download as pre-built binaries"
    echo "   ✅ Build completes in 10-20 minutes"
    echo "   ✅ No permission errors"
    echo ""
    echo "3. If builds still fail, check:"
    echo "   - Network connectivity to rstats-on-nix.cachix.org"
    echo "   - Disk space in /grphome/grp_batch_effects/nix/"
    echo ""
    echo "=========================================="
else
    echo ""
    echo "❌ Failed to configure cache with cachix"
    echo ""
    echo "This may require sudo privileges. Try:"
    echo "  sudo cachix use rstats-on-nix"
    echo ""
    echo "Or manually configure ~/.config/nix/nix.conf with:"
    echo "  substituters = https://cache.nixos.org https://rstats-on-nix.cachix.org"
    echo "  trusted-public-keys = cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= rstats-on-nix.cachix.org-1:dtPhHsUZNTNBdceD5/1SsWJ7p7Kv6hvzZvFzYmxl1yY="
    exit 1
fi
