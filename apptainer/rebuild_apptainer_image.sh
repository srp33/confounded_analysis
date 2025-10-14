#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# Load the Apptainer module (needed on supercomputer)
module load apptainer

# Build the Apptainer image directly on the host system
# (can't build Apptainer images from within Apptainer containers)
apptainer build --force combatseq_image.sif ~/confounded_analysis/apptainer/apptainer_combatseq.def

# Move the new image to overwrite the old one
SHARED_DIR=~/groups/grp_batch_effects
mv combatseq_image.sif "$SHARED_DIR/combatseq_image.sif"

echo "Successfully built and deployed new Apptainer image to $SHARED_DIR/combatseq_image.sif"