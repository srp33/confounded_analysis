#!/bin/bash

# BatchQC Analysis Wrapper Script
# Runs BatchQC on unadjusted and adjusted datasets

set -e

# Configuration
DATA_DIR="grp_batch_effects/data/paired_datasets"
OUTPUT_DIR="grp_batch_effects/outputs/batchqc"
SCRIPT_DIR="scripts/evaluations"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Function to run BatchQC on a single dataset
run_batchqc() {
    local dataset=$1
    local data_file=$2
    local suffix=$3
    
    echo "Running BatchQC on $dataset ($suffix)..."
    
    local output_subdir="$OUTPUT_DIR/${dataset}_${suffix}"
    
    if [[ -f "$data_file" ]]; then
        Rscript "$SCRIPT_DIR/batchqc_analysis.R" "$data_file" "$output_subdir"
        echo "✓ Completed: $dataset ($suffix)"
    else
        echo "⚠ File not found: $data_file"
    fi
}

# Find all datasets and run BatchQC
echo "Starting BatchQC analysis..."

for dataset_dir in "$DATA_DIR"/*/; do
    if [[ -d "$dataset_dir" ]]; then
        dataset=$(basename "$dataset_dir")
        echo "Processing dataset: $dataset"
        
        # Run on unadjusted data
        unadjusted_file="$dataset_dir/unadjusted.csv"
        run_batchqc "$dataset" "$unadjusted_file" "unadjusted"
        
        # Run on adjusted data (example adjusters)
        for adjuster in "gmm_affine" "gmm_nonlinear_unit_var" "combat" "limma"; do
            adjusted_file="$dataset_dir/${adjuster}.csv"
            if [[ -f "$adjusted_file" ]]; then
                run_batchqc "$dataset" "$adjusted_file" "$adjuster"
            fi
        done
        
        echo "---"
    fi
done

echo "BatchQC analysis complete!"
echo "Reports saved in: $OUTPUT_DIR"