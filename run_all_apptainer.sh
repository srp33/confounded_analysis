#! /bin/bash
./init_apptainer.sh
apptainer exec $APPTAINER_IMAGE /scripts/all.sh
