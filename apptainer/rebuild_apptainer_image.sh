#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# Load the Apptainer module (needed on supercomputer)
module load apptainer

# Build the Apptainer image directly on the host system
# (can't build Apptainer images from within Apptainer containers)
apptainer build --force remove-batch-effects.sif apptainer.def

# Move the new image to overwrite the old one
SHARED_DIR=~/groups/grp_batch_effects
mv remove-batch-effects.sif "$SHARED_DIR/remove-batch-effects.sif"

echo "Successfully built and deployed new Apptainer image to $SHARED_DIR/remove-batch-effects.sif"