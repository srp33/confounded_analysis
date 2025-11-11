#!/bin/bash
# run_generator.sh
# Phase 1 (Authoring): Run generate_env.R to create default.nix and .Rprofile
# This script uses a temporary R+rix shell to generate the environment definition

set -euo pipefail

# Configuration
NIX_ROOT="/grphome/grp_batch_effects/nix"
NIX_CHROOT_CMD="$NIX_ROOT/nix-user-chroot"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Phase 1: Authoring (Generating ComBat-seq Environment Definition) ===${NC}"
echo ""
echo "This script will:"
echo "  1. Enter a temporary Nix shell with R and rix"
echo "  2. Run generate_env.R to create default.nix and .Rprofile"
echo "  3. Exit the temporary shell"
echo ""
echo "Location: $SCRIPT_DIR"
echo "Nix root: $NIX_ROOT"
echo "Target: Bioconductor 3.11 (R 4.0.x era)"
echo ""

# Check if nix-user-chroot exists
if [ ! -f "$NIX_CHROOT_CMD" ]; then
    echo -e "${RED}ERROR: nix-user-chroot not found at $NIX_CHROOT_CMD${NC}"
    echo "Please ensure Nix is installed in the group directory."
    echo "See Task 1 in the migration plan for installation instructions."
    exit 1
fi

# Check if generate_env.R exists
if [ ! -f "$SCRIPT_DIR/generate_env.R" ]; then
    echo -e "${RED}ERROR: generate_env.R not found in $SCRIPT_DIR${NC}"
    exit 1
fi

# Change to script directory (important for rix to generate files in correct location)
cd "$SCRIPT_DIR"

echo -e "${YELLOW}Entering temporary authoring shell...${NC}"
echo "This may take a minute if R and rix need to be downloaded."
echo ""

# Run Phase 1: Enter nix-user-chroot, source nix profile, and run generate_env.R
# Note: We'll install rix from CRAN in the R session since rPackages.rix has build issues
$NIX_CHROOT_CMD "$NIX_ROOT" bash -c "
    # Source Nix profile to make nix-shell available
    source ~/.nix-profile/etc/profile.d/nix.sh
    
    echo 'Running generate_env.R in temporary R shell...'
    echo 'Installing rix package from CRAN if needed...'
    echo ''
    
    # Enter temporary shell with R and system dependencies, then run the helper script
    cd '$SCRIPT_DIR' && nix-shell -p R curl openssl libxml2 --run 'Rscript run_with_rix.R'
"

# Check if files were created
echo ""
echo -e "${GREEN}=== Checking generated files ===${NC}"

if [ -f "$SCRIPT_DIR/default.nix" ]; then
    echo -e "${GREEN}✓ default.nix created successfully${NC}"
    echo "  Size: $(du -h "$SCRIPT_DIR/default.nix" | cut -f1)"
    echo "  Lines: $(wc -l < "$SCRIPT_DIR/default.nix")"
else
    echo -e "${RED}✗ default.nix not found${NC}"
    exit 1
fi

if [ -f "$SCRIPT_DIR/.Rprofile" ]; then
    echo -e "${GREEN}✓ .Rprofile created successfully${NC}"
    echo "  Size: $(du -h "$SCRIPT_DIR/.Rprofile" | cut -f1)"
    echo "  Lines: $(wc -l < "$SCRIPT_DIR/.Rprofile")"
else
    echo -e "${RED}✗ .Rprofile not found${NC}"
    exit 1
fi

echo ""
echo -e "${GREEN}=== Phase 1 Complete ===${NC}"
echo ""
echo "ComBat-seq environment definition created successfully!"
echo ""
echo "Next steps:"
echo "  1. Review default.nix to verify package list:"
echo "     less default.nix"
echo ""
echo "  2. Build and test the environment (Phase 2):"
echo "     See Task 8.4 for Phase 2 instructions"
echo ""
echo "  3. To rebuild after modifying generate_env.R:"
echo "     ./run_generator.sh"
echo ""
