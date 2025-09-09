#! /bin/bash

set -e

SHARED_DIR=~/groups/grp_batch_effects
ANALYSIS_DIR=~/confounded_analysis
SCRIPTS_DIR=$ANALYSIS_DIR/scripts

# Make directories that need to exist
mkdir -p $SHARED_DIR/data $SHARED_DIR/data/.cache 
mkdir -p $SHARED_DIR/outputs/figures/ $SHARED_DIR/outputs/metrics $SHARED_DIR/outputs/tables

# Specific to supercomputer, loads the apptainer executable
module load apptainer 

# Export creates global variables that are useful later on
export APPTAINER_BINDPATH="$SHARED_DIR/data:/data,$SHARED_DIR/outputs:/outputs,$SCRIPTS_DIR:/scripts"
export APPTAINER_IMAGE="$SHARED_DIR/remove-batch-effects.sif"