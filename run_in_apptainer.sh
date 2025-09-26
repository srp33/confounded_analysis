#! /bin/bash

# Check if an argument was provided
if [ -z "$1" ]; then
    echo "Usage: $0 <executable>"
    exit 1
fi

# Assign the first argument to a variable
SCRIPT_TO_RUN=$1
shift  # remove the first arg, leaving only the rest

# Source the initialization script
source ~/confounded_analysis/init_apptainer.sh

# Execute the specified script inside the Apptainer container, forwarding all args
newgrp grp_batch_effects
apptainer exec --contain "$APPTAINER_IMAGE" "$SCRIPT_TO_RUN" "$@"