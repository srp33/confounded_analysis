#!/bin/bash
#SBATCH --job-name=prepare_omics
#SBATCH --output=logs/omics_%A_%a.log
#SBATCH --error=logs/omics_%A_%a.log
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G


DATA_DIR="/grphome/grp_batch_effects/outputs/metadata_features/labeled_adjusted"
SCRIPT=format_gseo_input.py
META_PREFIX=meta_

mkdir -p logs tmp_lists

# ---------------------------------------------------
# STEP 1: Generate class column list from first CSV
# ---------------------------------------------------

FIRST_FILE=$(find $DATA_DIR -name "*.csv" | head -n 1)

pixi run python - <<EOF > tmp_lists/class_columns.txt
import pandas as pd
df = pd.read_csv("$FIRST_FILE", nrows=1)
cols = [c for c in df.columns if c.startswith("$META_PREFIX")]
for c in cols:
    print(c)
EOF

# ---------------------------------------------------
# STEP 2: Generate input file list (one per subfolder)
# ---------------------------------------------------

find $DATA_DIR -mindepth 2 -maxdepth 2 -name "*.csv" > tmp_lists/input_files.txt

NUM_CLASSES=$(wc -l < tmp_lists/class_columns.txt)
NUM_FILES=$(wc -l < tmp_lists/input_files.txt)

TOTAL=$((NUM_CLASSES * NUM_FILES))

# ---------------------------------------------------
# STEP 3: Run array job
# ---------------------------------------------------

if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Submit with:"
    echo "sbatch --array=1-$TOTAL run_prepare_omics.slurm"
    exit 0
fi

TASK_ID=$((SLURM_ARRAY_TASK_ID - 1))

FILE_INDEX=$((TASK_ID / NUM_CLASSES + 1))
CLASS_INDEX=$((TASK_ID % NUM_CLASSES + 1))

INPUT_FILE=$(sed -n "${FILE_INDEX}p" tmp_lists/input_files.txt)
CLASS_COLUMN=$(sed -n "${CLASS_INDEX}p" tmp_lists/class_columns.txt)

DATASET=$(basename $(dirname $INPUT_FILE))
CLASS_NAME=$(echo $CLASS_COLUMN | sed 's/meta_//')

OUT_PREFIX=${DATASET}_${CLASS_NAME}

echo "Processing:"
echo "  Input: $INPUT_FILE"
echo "  Class column: $CLASS_COLUMN"
echo "  Output prefix: $OUT_PREFIX"

pixi run python $SCRIPT "$INPUT_FILE" \
    --class_column "$CLASS_COLUMN" \
    --sample_column meta_source \
    --out_prefix "$OUT_PREFIX"