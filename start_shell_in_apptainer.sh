#! /bin/bash
# Source the initialization script
source ~/confounded_analysis/init_apptainer.sh

# Execute the specified script inside the Apptainer container, forwarding all args
newgrp grp_batch_effects
apptainer shell "$APPTAINER_IMAGE"