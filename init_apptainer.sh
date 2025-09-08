#! /bin/bash

set -e

SHARED_DIR=~/groups/grp_batch_effects

mkdir -p $SHARED_DIR/data $SHARED_DIR/data/.cache 
mkdir -p $SHARED_DIR/outputs/figures/ $SHARED_DIR/outputs/metrics $SHARED_DIR/outputs/tables

module load apptainer

export APPTAINER_BINDPATH="$SHARED_DIR/data:/data,$SHARED_DIR/outputs:/outputs,$(pwd)/scripts:/scripts"
export APPTAINER_IMAGE=/groups/grp_batch_effects/outputs/images/batch_effects.sif