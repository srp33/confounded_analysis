#! /bin/bash

set -e

image="remove-batch-effects.sif"

SHARED_DIR=$(pwd)/../groups/grp_batch_effects

echo "Shared dir:"
echo $SHARED_DIR

mkdir -p $SHARED_DIR/data $SHARED_DIR/data/.cache 
mkdir -p $SHARED_DIR/outputs/figures/ $SHARED_DIR/outputs/metrics $SHARED_DIR/outputs/tables
mkdir -p $SHARED_DIR/tmp

module load apptainer

apptainer exec \
  --bind "$SHARED_DIR/data":/data \
  --bind "$SHARED_DIR/outputs":/outputs \
  --bind "$(pwd)/scripts":/scripts \
  --bind "$SHARED_DIR/data/.cache":/.cache \
  $SHARED_DIR/$image \
  /scripts/all.sh
