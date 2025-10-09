#!/bin/bash

# Complete GMM Integration Test Script
# This script runs a comprehensive test of the GMM integration in the robustifying analysis

echo "=========================================="
echo "Complete GMM Integration Test"
echo "Date: $(date)"
echo "=========================================="

# Set working directory
cd /scripts/evaluations/robustifying

# Test parameters for comprehensive validation
TEST_PARAMS=(
    "20 0 1"  # Basic test case
    "20 3 2"  # Medium batch effect
    "40 0 1"  # Larger sample size
)

echo ""
echo "Test parameters to validate:"
for params in "${TEST_PARAMS[@]}"; do
    echo "  - N_sample_size=$(echo $params | cut -d' ' -f1), max_batch_mean=$(echo $params | cut -d' ' -f2), max_batch_var=$(echo $params | cut -d' ' -f3)"
done

echo ""
echo "=========================================="
echo "1. Validating Prerequisites"
echo "=========================================="

# Check if gmm_adjust.R exists
if [ -f "../../adjust/gmm_adjust.R" ]; then
    echo "✓ gmm_adjust.R found"
else
    echo "✗ gmm_adjust.R not found"
    exit 1
fi

# Check if data file exists
if [ -f "data/combined_sub.RData" ]; then
    echo "✓ Data file found"
else
    echo "✗ Data file not found"
    exit 1
fi

# Check if 1_simpipe.R exists and has GMM integration
if grep -q "gmm_adjust" code/1_simpipe.R; then
    echo "✓ GMM integration found in 1_simpipe.R"
else
    echo "✗ GMM integration not found in 1_simpipe.R"
    exit 1
fi

echo ""
echo "=========================================="
echo "2. Running Simulation Tests"
echo "=========================================="

SUCCESS_COUNT=0
TOTAL_TESTS=${#TEST_PARAMS[@]}

for i in "${!TEST_PARAMS[@]}"; do
    params=${TEST_PARAMS[$i]}
    N_SAMPLE=$(echo $params | cut -d' ' -f1)
    M_BATCH=$(echo $params | cut -d' ' -f2)
    V_BATCH=$(echo $params | cut -d' ' -f3)
    
    echo ""
    echo "Test $((i+1))/$TOTAL_TESTS: Running simulation with parameters $params"
    echo "----------------------------------------"
    
    # Run the simulation
    echo "Executing: Rscript code/1_simpipe.R $N_SAMPLE $M_BATCH $V_BATCH"
    
    if timeout 300 Rscript code/1_simpipe.R $N_SAMPLE $M_BATCH $V_BATCH; then
        echo "✓ Simulation completed successfully"
        
        # Check output files for each learner type
        LEARNERS=("lasso" "rf" "svm")
        METRICS=("auc" "mxe")
        
        for learner in "${LEARNERS[@]}"; do
            for metric in "${METRICS[@]}"; do
                OUTPUT_FILE="results/${learner}_${metric}_batchN${N_SAMPLE}_m${M_BATCH}_v${V_BATCH}.csv"
                
                if [ -f "$OUTPUT_FILE" ]; then
                    echo "  ✓ Output file created: $OUTPUT_FILE"
                    
                    # Check if GMM column exists
                    if head -1 "$OUTPUT_FILE" | grep -q "GMM"; then
                        echo "    ✓ GMM column found"
                        
                        # Check if GMM has actual values (not just header)
                        if [ $(wc -l < "$OUTPUT_FILE") -gt 1 ]; then
                            GMM_VALUES=$(tail -n +2 "$OUTPUT_FILE" | cut -d',' -f7 | grep -v "^$" | wc -l)
                            if [ $GMM_VALUES -gt 0 ]; then
                                echo "    ✓ GMM has $GMM_VALUES data points"
                            else
                                echo "    ✗ GMM column exists but has no data"
                            fi
                        fi
                    else
                        echo "    ✗ GMM column not found in $OUTPUT_FILE"
                        echo "    Available columns: $(head -1 "$OUTPUT_FILE")"
                    fi
                else
                    echo "  ✗ Output file not created: $OUTPUT_FILE"
                fi
            done
        done
        
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "✗ Simulation failed or timed out"
    fi
done

echo ""
echo "=========================================="
echo "3. Validation Summary"
echo "=========================================="

echo "Successful tests: $SUCCESS_COUNT/$TOTAL_TESTS"

if [ $SUCCESS_COUNT -eq $TOTAL_TESTS ]; then
    echo "✓ All tests passed - GMM integration is working correctly"
    
    echo ""
    echo "Sample output from latest test:"
    LATEST_FILE=$(ls -t results/lasso_auc_*.csv | head -1)
    if [ -f "$LATEST_FILE" ]; then
        echo "File: $LATEST_FILE"
        echo "Header: $(head -1 "$LATEST_FILE")"
        echo "Sample data: $(tail -1 "$LATEST_FILE")"
    fi
    
    exit 0
else
    echo "✗ Some tests failed - GMM integration needs attention"
    exit 1
fi

echo ""
echo "=========================================="
echo "Test Complete"
echo "=========================================="