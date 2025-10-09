#!/bin/bash

# Test script for GMM integration in robustifying analysis
# This script runs a single simulation test to validate the GMM integration

echo "=== GMM Integration Test ==="
echo "Testing GMM integration in robustifying analysis pipeline"
echo "Date: $(date)"
echo ""

# Set test parameters (small values for quick testing)
N_SAMPLE_SIZE=20
MAX_BATCH_MEAN=0
MAX_BATCH_VAR=1

echo "Test parameters:"
echo "  N_sample_size: $N_SAMPLE_SIZE"
echo "  max_batch_mean: $MAX_BATCH_MEAN"
echo "  max_batch_var: $MAX_BATCH_VAR"
echo ""

# Change to robustifying directory
cd /scripts/evaluations/robustifying

echo "Current directory: $(pwd)"
echo "Available files:"
ls -la code/
echo ""

# Check if gmm_adjust.R is available
if [ -f "../../adjust/gmm_adjust.R" ]; then
    echo "✓ gmm_adjust.R found at ../../adjust/gmm_adjust.R"
else
    echo "✗ gmm_adjust.R not found at ../../adjust/gmm_adjust.R"
    exit 1
fi

# Check if data file exists
if [ -f "data/combined_sub.RData" ]; then
    echo "✓ Data file found: data/combined_sub.RData"
else
    echo "✗ Data file not found: data/combined_sub.RData"
    exit 1
fi

echo ""
echo "=== Running simulation test ==="

# Run the simulation with test parameters
Rscript code/1_simpipe.R $N_SAMPLE_SIZE $MAX_BATCH_MEAN $MAX_BATCH_VAR

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "=== Test completed successfully ==="
    
    # Check if output files were created
    echo "Checking output files:"
    
    # Expected output file pattern
    OUTPUT_FILE="results/lasso_auc_batchN${N_SAMPLE_SIZE}_m${MAX_BATCH_MEAN}_v${MAX_BATCH_VAR}.csv"
    if [ -f "$OUTPUT_FILE" ]; then
        echo "✓ Output file created: $OUTPUT_FILE"
        echo "File contents:"
        head -5 "$OUTPUT_FILE"
        
        # Check if GMM column exists
        if grep -q "GMM" "$OUTPUT_FILE"; then
            echo "✓ GMM column found in output"
        else
            echo "✗ GMM column not found in output"
        fi
    else
        echo "✗ Expected output file not found: $OUTPUT_FILE"
        echo "Available files in results/:"
        ls -la results/ | head -10
    fi
    
else
    echo ""
    echo "=== Test failed ==="
    echo "Exit code: $?"
    exit 1
fi

echo ""
echo "=== GMM Integration Test Complete ==="