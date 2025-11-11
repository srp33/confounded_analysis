#!/bin/bash
# test_task_8.3.sh
# Comprehensive test for Task 8.3: Build and test batch-effects environment

set -euo pipefail

NIX_ROOT="/grphome/grp_batch_effects/nix"
ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_FILE="$ENV_DIR/test_task_8.3.log"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "$1" | tee -a "$LOG_FILE"
}

# Test result tracking
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_TOTAL=0

# Test function
run_test() {
    local test_name="$1"
    local test_command="$2"
    
    TESTS_TOTAL=$((TESTS_TOTAL + 1))
    log "\n${YELLOW}Test $TESTS_TOTAL: $test_name${NC}"
    
    if eval "$test_command" >> "$LOG_FILE" 2>&1; then
        log "${GREEN}✓ PASSED${NC}"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        log "${RED}✗ FAILED${NC}"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# Initialize log
echo "Task 8.3 Test Results - $(date)" > "$LOG_FILE"
echo "======================================" >> "$LOG_FILE"

log "${YELLOW}=== Task 8.3: Build and Test batch-effects Environment ===${NC}\n"

# Test 1: Check R version
run_test "Check R version" \
    "$NIX_ROOT/nix-user-chroot $NIX_ROOT bash -c \"
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        cd $ENV_DIR && \
        timeout 600 nix-shell --option substituters 'https://cache.nixos.org https://rstats-on-nix.cachix.org' --option trusted-public-keys 'cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= rstats-on-nix.cachix.org-1:vdiiVgocg6WeJrODIqdprZRUrhi1JzhBnXv7aWI6+F0=' --run 'R --version | head -1'
    \""

# Test 2: Test tidyverse loading
run_test "Test tidyverse package loading" \
    "$NIX_ROOT/nix-user-chroot $NIX_ROOT bash -c \"
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        cd $ENV_DIR && \
        timeout 600 nix-shell --option substituters 'https://cache.nixos.org https://rstats-on-nix.cachix.org' --option trusted-public-keys 'cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= rstats-on-nix.cachix.org-1:vdiiVgocg6WeJrODIqdprZRUrhi1JzhBnXv7aWI6+F0=' --run 'Rscript -e \\\"library(tidyverse); cat(\\\\\\\"tidyverse loaded successfully\\\\n\\\\\\\")\\\"'
    \""

# Test 3: Verify .Rprofile isolation
run_test "Verify .Rprofile isolation (no system R contamination)" \
    "$NIX_ROOT/nix-user-chroot $NIX_ROOT bash -c \"
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        cd $ENV_DIR && \
        timeout 600 nix-shell --option substituters 'https://cache.nixos.org https://rstats-on-nix.cachix.org' --option trusted-public-keys 'cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= rstats-on-nix.cachix.org-1:vdiiVgocg6WeJrODIqdprZRUrhi1JzhBnXv7aWI6+F0=' --run 'Rscript -e \\\".libPaths()\\\" | grep -q /nix/store'
    \""

# Test 4: Test sample R script
log "\n${YELLOW}Test $((TESTS_TOTAL + 1)): Test with sample R script${NC}"
TESTS_TOTAL=$((TESTS_TOTAL + 1))

cat > "$ENV_DIR/test_script.R" << 'EOF'
# Simple test script
library(dplyr)
library(ggplot2)

# Create test data
df <- data.frame(
    x = 1:10,
    y = rnorm(10)
)

# Test dplyr
result <- df %>%
    filter(x > 5) %>%
    summarise(mean_y = mean(y))

cat("Test script completed successfully\n")
cat("Mean of filtered data:", result$mean_y, "\n")
EOF

if $NIX_ROOT/nix-user-chroot $NIX_ROOT bash -c "
    source ~/.nix-profile/etc/profile.d/nix.sh && \
    cd $ENV_DIR && \
    timeout 600 nix-shell --option substituters 'https://cache.nixos.org https://rstats-on-nix.cachix.org' --option trusted-public-keys 'cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= rstats-on-nix.cachix.org-1:vdiiVgocg6WeJrODIqdprZRUrhi1JzhBnXv7aWI6+F0=' --run 'Rscript test_script.R'
" >> "$LOG_FILE" 2>&1; then
    log "${GREEN}✓ PASSED${NC}"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    log "${RED}✗ FAILED${NC}"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

# Test 5: Measure Nix store size
log "\n${YELLOW}Test $((TESTS_TOTAL + 1)): Measure Nix store size${NC}"
TESTS_TOTAL=$((TESTS_TOTAL + 1))

if STORE_SIZE=$(du -sh /grphome/grp_batch_effects/nix/store/ 2>/dev/null | cut -f1); then
    log "${GREEN}✓ PASSED${NC}"
    log "Nix store size: $STORE_SIZE"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    log "${RED}✗ FAILED${NC}"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

# Summary
log "\n${YELLOW}=== Test Summary ===${NC}"
log "Total tests: $TESTS_TOTAL"
log "${GREEN}Passed: $TESTS_PASSED${NC}"
log "${RED}Failed: $TESTS_FAILED${NC}"

if [ $TESTS_FAILED -eq 0 ]; then
    log "\n${GREEN}✓ All tests passed! Task 8.3 is complete.${NC}"
    exit 0
else
    log "\n${RED}✗ Some tests failed. See log for details: $LOG_FILE${NC}"
    exit 1
fi
