#!/bin/bash
# test_environment.sh
# Test script for Phase 2: Build and validate the batch-effects R environment

set -euo pipefail

# Configuration
NIX_ROOT="/grphome/grp_batch_effects/nix"
NIX_CHROOT_CMD="$NIX_ROOT/nix-user-chroot"
ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=========================================="
echo "Testing batch-effects R Environment"
echo "=========================================="
echo ""
echo "Environment directory: $ENV_DIR"
echo "Nix root: $NIX_ROOT"
echo ""

# Check prerequisites
echo "=== Checking prerequisites ==="
if [ ! -f "$NIX_CHROOT_CMD" ]; then
    echo "ERROR: nix-user-chroot not found at $NIX_CHROOT_CMD"
    exit 1
fi
echo "✓ nix-user-chroot found"

if [ ! -f "$ENV_DIR/default.nix" ]; then
    echo "ERROR: default.nix not found"
    exit 1
fi
echo "✓ default.nix found"

if [ ! -f "$ENV_DIR/.Rprofile" ]; then
    echo "WARNING: .Rprofile not found"
else
    echo "✓ .Rprofile found"
fi
echo ""

# Function to run commands in nix-shell
run_in_nix() {
    local cmd="$1"
    cd "$ENV_DIR"
    "$NIX_CHROOT_CMD" "$NIX_ROOT" bash -c "
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        nix-shell --run '$cmd'
    "
}

# Test 1: Check R version
echo "=== Test 1: Check R version ==="
echo "Running: nix-shell --run 'R --version'"
START_TIME=$(date +%s)
if run_in_nix "R --version" 2>&1 | head -5; then
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    echo "✓ R version check passed (${DURATION}s)"
else
    echo "✗ R version check failed"
    exit 1
fi
echo ""

# Test 2: Test tidyverse loading
echo "=== Test 2: Test tidyverse package loading ==="
echo "Running: nix-shell --run 'Rscript -e \"library(tidyverse)\"'"
if run_in_nix "Rscript -e 'library(tidyverse)'" 2>&1 | tail -10; then
    echo "✓ tidyverse loaded successfully"
else
    echo "✗ tidyverse loading failed"
    exit 1
fi
echo ""

# Test 3: Test Bioconductor packages
echo "=== Test 3: Test Bioconductor packages ==="
echo "Running: nix-shell --run 'Rscript -e \"library(limma); library(sva)\"'"
if run_in_nix "Rscript -e 'library(limma); library(sva); cat(\"Bioconductor packages loaded\\n\")'" 2>&1 | tail -5; then
    echo "✓ Bioconductor packages loaded successfully"
else
    echo "✗ Bioconductor package loading failed"
    exit 1
fi
echo ""

# Test 4: Test machine learning packages
echo "=== Test 4: Test machine learning packages ==="
echo "Running: nix-shell --run 'Rscript -e \"library(caret); library(xgboost)\"'"
if run_in_nix "Rscript -e 'library(caret); library(xgboost); cat(\"ML packages loaded\\n\")'" 2>&1 | tail -5; then
    echo "✓ Machine learning packages loaded successfully"
else
    echo "✗ Machine learning package loading failed"
    exit 1
fi
echo ""

# Test 5: Verify .Rprofile isolation
echo "=== Test 5: Verify .Rprofile isolation ==="
echo "Checking that system R libraries are not in library path..."
if run_in_nix "Rscript -e 'paths <- .libPaths(); cat(\"Library paths:\\n\"); print(paths); if(any(grepl(\"R_LIBS_USER\", paths))) stop(\"System R contamination detected\") else cat(\"✓ No system R contamination\\n\")'" 2>&1 | tail -10; then
    echo "✓ Library path isolation verified"
else
    echo "✗ Library path isolation check failed"
    exit 1
fi
echo ""

# Test 6: Test with a simple R script
echo "=== Test 6: Test with sample R script ==="
cat > /tmp/test_batch_effects.R << 'RSCRIPT'
# Simple test script to verify environment functionality
library(tidyverse)
library(limma)
library(sva)

# Create test data
set.seed(42)
n_genes <- 100
n_samples <- 20
expr_data <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
batch <- rep(1:2, each = n_samples/2)

# Test ComBat
cat("Testing ComBat batch correction...\n")
corrected <- ComBat(dat = expr_data, batch = batch, mod = NULL)
cat("✓ ComBat completed successfully\n")

# Test basic tidyverse operations
cat("Testing tidyverse operations...\n")
df <- tibble(x = 1:10, y = x^2)
result <- df %>% mutate(z = x + y) %>% summarize(mean_z = mean(z))
cat("✓ Tidyverse operations completed\n")

cat("\n=== All R functionality tests passed ===\n")
RSCRIPT

echo "Running sample R script..."
if run_in_nix "Rscript /tmp/test_batch_effects.R" 2>&1 | tail -15; then
    echo "✓ Sample R script executed successfully"
else
    echo "✗ Sample R script execution failed"
    exit 1
fi
rm -f /tmp/test_batch_effects.R
echo ""

# Measure Nix store size
echo "=== Measuring Nix store size ==="
if [ -d "$NIX_ROOT/store" ]; then
    STORE_SIZE=$(du -sh "$NIX_ROOT/store" 2>/dev/null | cut -f1 || echo "unknown")
    echo "Nix store size: $STORE_SIZE"
else
    echo "Nix store not found at $NIX_ROOT/store"
fi
echo ""

# Summary
echo "=========================================="
echo "✓ All tests passed!"
echo "=========================================="
echo ""
echo "The batch-effects R environment is ready to use."
echo ""
echo "Usage examples:"
echo "  1. Interactive shell:"
echo "     cd $ENV_DIR && ./nix_activate.sh"
echo ""
echo "  2. Run R script:"
echo "     cd $ENV_DIR && $NIX_CHROOT_CMD $NIX_ROOT bash -c \\"
echo "       'source ~/.nix-profile/etc/profile.d/nix.sh && \\"
echo "        nix-shell --run \"Rscript /path/to/script.R\"'"
echo ""
echo "  3. Check installed packages:"
echo "     cd $ENV_DIR && $NIX_CHROOT_CMD $NIX_ROOT bash -c \\"
echo "       'source ~/.nix-profile/etc/profile.d/nix.sh && \\"
echo "        nix-shell --run \"Rscript -e \\\"installed.packages()[,c(\\\\\\\"Package\\\\\\\",\\\\\\\"Version\\\\\\\")]\\\"\"'"
echo ""
