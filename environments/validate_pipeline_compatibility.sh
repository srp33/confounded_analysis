#!/bin/bash
# Pipeline Script Compatibility Validation
# Tests representative pipeline scripts to ensure they work with uv/rix environments
#
# Usage:
#   bash environments/validate_pipeline_compatibility.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(dirname "$SCRIPT_DIR")"
SCRIPTS_DIR="${ANALYSIS_DIR}/scripts"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test results
TESTS_PASSED=0
TESTS_FAILED=0

print_header() {
    echo ""
    echo "======================================================================"
    echo "  $1"
    echo "======================================================================"
}

print_section() {
    echo ""
    echo "--- $1 ---"
}

print_success() {
    echo -e "${GREEN}✓${NC} $1"
    TESTS_PASSED=$((TESTS_PASSED + 1))
}

print_failure() {
    echo -e "${RED}✗${NC} $1"
    TESTS_FAILED=$((TESTS_FAILED + 1))
}

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

# ============================================================================
# Test 1: Python Data Processing Script
# ============================================================================
test_python_data_processing() {
    print_section "Test 1: Python Data Processing Script"
    
    # Find a representative Python script
    local test_scripts=(
        "${SCRIPTS_DIR}/prepdata/combine_datasets.py"
        "${SCRIPTS_DIR}/prepdata/preview_data.py"
        "${SCRIPTS_DIR}/evaluations/metrics.py"
    )
    
    local test_script=""
    for script in "${test_scripts[@]}"; do
        if [ -f "$script" ]; then
            test_script="$script"
            break
        fi
    done
    
    if [ -z "$test_script" ]; then
        print_warning "No suitable Python data processing script found, skipping"
        return 0
    fi
    
    echo "  Testing script: $test_script"
    
    # Test with --help flag (should not require data files)
    if ${SCRIPT_DIR}/run_with_env.sh "$test_script" --help &>/dev/null; then
        print_success "Python script help flag works"
    else
        # Some scripts might not have --help, try -h
        if ${SCRIPT_DIR}/run_with_env.sh "$test_script" -h &>/dev/null; then
            print_success "Python script help flag works"
        else
            print_warning "Python script does not support --help flag (not a failure)"
        fi
    fi
    
    # Test import validation (create a simple test)
    local import_test=$(mktemp /tmp/test_imports.XXXXXX.py)
    cat > "$import_test" << 'EOF'
#!/usr/bin/env python3
# Test that pipeline modules can be imported
import sys
from pathlib import Path

# Add scripts directory to path
scripts_dir = Path(__file__).parent.parent / "scripts"
sys.path.insert(0, str(scripts_dir))

# Test imports
try:
    from utils import DataFrameCache
    print("✓ Imported utils.DataFrameCache")
    
    from evaluations.util import repeated_cross_val
    print("✓ Imported evaluations.util.repeated_cross_val")
    
    from prepdata.config import DATASETS_DIR
    print("✓ Imported prepdata.config.DATASETS_DIR")
    
    print("All pipeline imports successful")
except ImportError as e:
    print(f"Import failed: {e}")
    sys.exit(1)
EOF
    
    if ${SCRIPT_DIR}/run_with_env.sh python "$import_test"; then
        print_success "Pipeline module imports work correctly"
    else
        print_failure "Pipeline module imports failed"
    fi
    
    rm -f "$import_test"
}

# ============================================================================
# Test 2: R Analysis Script
# ============================================================================
test_r_analysis_script() {
    print_section "Test 2: R Analysis Script"
    
    # Find a representative R script
    local test_scripts=(
        "${SCRIPTS_DIR}/adjust/gmm_adjust.R"
        "${SCRIPTS_DIR}/adjust/adjust.R"
        "${SCRIPTS_DIR}/evaluations/batchqc/batchqc_analysis.R"
    )
    
    local test_script=""
    for script in "${test_scripts[@]}"; do
        if [ -f "$script" ]; then
            test_script="$script"
            break
        fi
    done
    
    if [ -z "$test_script" ]; then
        print_warning "No suitable R analysis script found, skipping"
        return 0
    fi
    
    echo "  Testing script: $test_script"
    
    # Test with --help flag (should not require data files)
    if ${SCRIPT_DIR}/run_with_env.sh "$test_script" --help &>/dev/null; then
        print_success "R script help flag works"
    else
        # Some scripts might not have --help
        print_warning "R script does not support --help flag (not a failure)"
    fi
    
    # Test library loading (create a simple test)
    local library_test=$(mktemp /tmp/test_libraries.XXXXXX.R)
    cat > "$library_test" << 'EOF'
#!/usr/bin/env Rscript
# Test that pipeline libraries can be loaded

cat("Testing R library loading...\n")

# Core libraries used in pipeline
libraries <- c(
  "dplyr",
  "ggplot2",
  "tidyr",
  "glmnet",
  "caret",
  "sva",
  "limma"
)

failed <- character(0)

for (lib in libraries) {
  result <- tryCatch({
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
    cat(sprintf("✓ Loaded %s\n", lib))
    TRUE
  }, error = function(e) {
    cat(sprintf("✗ Failed to load %s: %s\n", lib, e$message))
    failed <<- c(failed, lib)
    FALSE
  })
}

if (length(failed) == 0) {
  cat("All pipeline libraries loaded successfully\n")
  quit(status = 0)
} else {
  cat(sprintf("Failed to load %d libraries: %s\n", length(failed), paste(failed, collapse = ", ")))
  quit(status = 1)
}
EOF
    
    if ${SCRIPT_DIR}/run_with_env.sh Rscript "$library_test"; then
        print_success "Pipeline R libraries load correctly"
    else
        print_failure "Pipeline R libraries failed to load"
    fi
    
    rm -f "$library_test"
}

# ============================================================================
# Test 3: Shell Script Wrapper
# ============================================================================
test_shell_script_wrapper() {
    print_section "Test 3: Shell Script Wrapper"
    
    # Find a representative shell script
    local test_scripts=(
        "${SCRIPTS_DIR}/adjust/all.sh"
        "${SCRIPTS_DIR}/evaluations/all.sh"
    )
    
    local test_script=""
    for script in "${test_scripts[@]}"; do
        if [ -f "$script" ]; then
            test_script="$script"
            break
        fi
    done
    
    if [ -z "$test_script" ]; then
        print_warning "No suitable shell script found, skipping"
        return 0
    fi
    
    echo "  Testing script: $test_script"
    
    # Create a simple test shell script
    local shell_test=$(mktemp /tmp/test_shell.XXXXXX.sh)
    cat > "$shell_test" << 'EOF'
#!/bin/bash
set -euo pipefail

echo "Testing shell script environment..."

# Check Python is available
if command -v python &>/dev/null; then
    echo "✓ Python available: $(python --version)"
else
    echo "✗ Python not available"
    exit 1
fi

# Check R is available
if command -v R &>/dev/null; then
    echo "✓ R available: $(R --version | head -1)"
else
    echo "✗ R not available"
    exit 1
fi

# Test Python import
python -c "import numpy; print('✓ NumPy import successful')"

# Test R library
Rscript -e "library(dplyr); cat('✓ dplyr load successful\n')"

echo "Shell script environment test completed"
EOF
    
    chmod +x "$shell_test"
    
    if ${SCRIPT_DIR}/run_with_env.sh "$shell_test"; then
        print_success "Shell script wrapper works correctly"
    else
        print_failure "Shell script wrapper failed"
    fi
    
    rm -f "$shell_test"
}

# ============================================================================
# Test 4: Command-Line Argument Handling
# ============================================================================
test_argument_handling() {
    print_section "Test 4: Command-Line Argument Handling"
    
    # Create test scripts with various argument patterns
    
    # Python script with argparse
    local py_test=$(mktemp /tmp/test_args.XXXXXX.py)
    cat > "$py_test" << 'EOF'
#!/usr/bin/env python3
import argparse

parser = argparse.ArgumentParser(description='Test argument handling')
parser.add_argument('--input', required=True, help='Input file')
parser.add_argument('--output', required=True, help='Output file')
parser.add_argument('--verbose', action='store_true', help='Verbose output')

args = parser.parse_args()

print(f"Input: {args.input}")
print(f"Output: {args.output}")
print(f"Verbose: {args.verbose}")
print("✓ Python argument handling successful")
EOF
    
    if ${SCRIPT_DIR}/run_with_env.sh python "$py_test" --input test.csv --output result.csv --verbose; then
        print_success "Python argparse handling works"
    else
        print_failure "Python argparse handling failed"
    fi
    
    rm -f "$py_test"
    
    # R script with argparse
    local r_test=$(mktemp /tmp/test_args.XXXXXX.R)
    cat > "$r_test" << 'EOF'
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description='Test argument handling')
parser$add_argument('--input', required=TRUE, help='Input file')
parser$add_argument('--output', required=TRUE, help='Output file')
parser$add_argument('--verbose', action='store_true', help='Verbose output')

args <- parser$parse_args()

cat(sprintf("Input: %s\n", args$input))
cat(sprintf("Output: %s\n", args$output))
cat(sprintf("Verbose: %s\n", args$verbose))
cat("✓ R argument handling successful\n")
EOF
    
    if ${SCRIPT_DIR}/run_with_env.sh Rscript "$r_test" --input test.csv --output result.csv --verbose; then
        print_success "R argparse handling works"
    else
        print_failure "R argparse handling failed"
    fi
    
    rm -f "$r_test"
}

# ============================================================================
# Test 5: File I/O Operations
# ============================================================================
test_file_io() {
    print_section "Test 5: File I/O Operations"
    
    # Create temporary test directory
    local test_dir=$(mktemp -d /tmp/test_io.XXXXXX)
    
    # Python file I/O test
    local py_io_test=$(mktemp /tmp/test_io.XXXXXX.py)
    cat > "$py_io_test" << EOF
#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pathlib import Path

test_dir = Path("$test_dir")

# Test CSV I/O
df = pd.DataFrame({
    'A': np.random.randn(100),
    'B': np.random.randn(100),
    'C': np.random.choice(['X', 'Y', 'Z'], 100)
})

csv_file = test_dir / "test.csv"
df.to_csv(csv_file, index=False)
df_loaded = pd.read_csv(csv_file)
assert df_loaded.shape == df.shape
print(f"✓ Python CSV I/O: {csv_file}")

# Test HDF5 I/O
h5_file = test_dir / "test.h5"
df.to_hdf(h5_file, key='data', mode='w')
df_loaded = pd.read_hdf(h5_file, key='data')
assert df_loaded.shape == df.shape
print(f"✓ Python HDF5 I/O: {h5_file}")

print("Python file I/O test completed")
EOF
    
    if ${SCRIPT_DIR}/run_with_env.sh python "$py_io_test"; then
        print_success "Python file I/O works correctly"
    else
        print_failure "Python file I/O failed"
    fi
    
    rm -f "$py_io_test"
    
    # R file I/O test
    local r_io_test=$(mktemp /tmp/test_io.XXXXXX.R)
    cat > "$r_io_test" << EOF
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(readr))

test_dir <- "$test_dir"

# Test CSV I/O
df <- data.frame(
  x = rnorm(100),
  y = rnorm(100),
  group = sample(c("A", "B", "C"), 100, replace = TRUE)
)

csv_file <- file.path(test_dir, "test_r.csv")
write_csv(df, csv_file)
df_loaded <- read_csv(csv_file, show_col_types = FALSE)
stopifnot(nrow(df_loaded) == nrow(df))
cat(sprintf("✓ R CSV I/O: %s\n", csv_file))

# Test RDS I/O
rds_file <- file.path(test_dir, "test.rds")
saveRDS(df, rds_file)
df_loaded <- readRDS(rds_file)
stopifnot(nrow(df_loaded) == nrow(df))
cat(sprintf("✓ R RDS I/O: %s\n", rds_file))

cat("R file I/O test completed\n")
EOF
    
    if ${SCRIPT_DIR}/run_with_env.sh Rscript "$r_io_test"; then
        print_success "R file I/O works correctly"
    else
        print_failure "R file I/O failed"
    fi
    
    rm -f "$r_io_test"
    
    # Clean up test directory
    rm -rf "$test_dir"
}

# ============================================================================
# Test 6: Actual Pipeline Script (if data available)
# ============================================================================
test_actual_pipeline_script() {
    print_section "Test 6: Actual Pipeline Script (Dry Run)"
    
    # This test is optional and only runs if test data is available
    print_warning "Actual pipeline script testing requires real data"
    print_warning "Skipping full pipeline test (would require data setup)"
    
    # Instead, we'll test that scripts can at least be parsed/validated
    local scripts_to_check=(
        "${SCRIPTS_DIR}/adjust/gmm_adjust.R"
        "${SCRIPTS_DIR}/prepdata/combine_datasets.py"
        "${SCRIPTS_DIR}/evaluations/metrics.py"
    )
    
    local syntax_errors=0
    
    for script in "${scripts_to_check[@]}"; do
        if [ ! -f "$script" ]; then
            continue
        fi
        
        echo "  Checking syntax: $(basename $script)"
        
        case "$script" in
            *.py)
                if python -m py_compile "$script" 2>/dev/null; then
                    echo "    ✓ Python syntax valid"
                else
                    echo "    ✗ Python syntax error"
                    syntax_errors=$((syntax_errors + 1))
                fi
                ;;
            *.R)
                if Rscript -e "parse('$script')" &>/dev/null; then
                    echo "    ✓ R syntax valid"
                else
                    echo "    ✗ R syntax error"
                    syntax_errors=$((syntax_errors + 1))
                fi
                ;;
        esac
    done
    
    if [ $syntax_errors -eq 0 ]; then
        print_success "Pipeline script syntax validation passed"
    else
        print_failure "Found $syntax_errors syntax errors in pipeline scripts"
    fi
}

# ============================================================================
# Main Function
# ============================================================================
main() {
    print_header "Pipeline Script Compatibility Validation"
    echo "  Started: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "  Scripts directory: $SCRIPTS_DIR"
    echo ""
    
    # Check prerequisites
    if [ ! -f "${SCRIPT_DIR}/run_with_env.sh" ]; then
        echo "Error: run_with_env.sh not found at ${SCRIPT_DIR}/run_with_env.sh"
        exit 1
    fi
    
    if [ ! -d "$SCRIPTS_DIR" ]; then
        echo "Error: Scripts directory not found at $SCRIPTS_DIR"
        exit 1
    fi
    
    # Run tests
    test_python_data_processing || true
    test_r_analysis_script || true
    test_shell_script_wrapper || true
    test_argument_handling || true
    test_file_io || true
    test_actual_pipeline_script || true
    
    # Summary
    print_header "Test Summary"
    
    local total_tests=$((TESTS_PASSED + TESTS_FAILED))
    echo ""
    echo "  Total tests:   $total_tests"
    echo "  Passed:        $TESTS_PASSED"
    echo "  Failed:        $TESTS_FAILED"
    echo ""
    
    if [ $TESTS_FAILED -eq 0 ]; then
        echo -e "${GREEN}✓ ALL PIPELINE COMPATIBILITY TESTS PASSED${NC}"
        echo "======================================================================"
        return 0
    else
        echo -e "${RED}✗ SOME PIPELINE COMPATIBILITY TESTS FAILED${NC}"
        echo "======================================================================"
        return 1
    fi
}

# Run main function
main
exit $?
