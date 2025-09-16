#! /bin/bash
# Source the initialization script
source ~/confounded_analysis/init_apptainer.sh

# Execute the specified script inside the Apptainer container, forwarding all args
apptainer shell "$APPTAINER_IMAGE"