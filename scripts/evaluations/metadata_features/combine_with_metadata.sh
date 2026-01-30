#!/bin/bash
#SBATCH --job-name=combine_gold      # Job name
#SBATCH --output=logs/combine/combine_gold_%j.log # Standard output + error log
#SBATCH --time=02:00:00              # Max run time (HH:MM:SS)
#SBATCH --mem=32G                     # Memory per node
#SBATCH --cpus-per-task=4            # Number of CPU cores

# Set input and output paths
INPUT_DIR="/grphome/grp_batch_effects/data/gold"
OUTPUT_FILE="/grphome/grp_batch_effects/outputs/metadata_features/combined/combined_labeled.csv"
COMBINE_SCRIPT="../classify_er_all_datasets/combine_all.py"

mkdir -p "$(dirname "$OUTPUT_FILE")"

# Run Python script
pixi run python "$COMBINE_SCRIPT" --input-dir "$INPUT_DIR" --output-file "$OUTPUT_FILE"
