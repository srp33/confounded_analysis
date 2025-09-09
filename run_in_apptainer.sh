#! /bin/bash

# Check if an argument was provided
if [ -z "$1" ]; then
    echo "Usage: $0 <bash_script_path.sh>"
    exit 1
fi

# Assign the first argument to a variable
SCRIPT_TO_RUN=$1

# Source the initialization script
source ~/confounded_analysis/init_apptainer.sh

# Execute the specified script inside the Apptainer container
apptainer exec "$APPTAINER_IMAGE" bash "$SCRIPT_TO_RUN"