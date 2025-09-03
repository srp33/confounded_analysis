#!/bin/bash

set -e

# Set PYTHONPATH for all Python scripts
export PYTHONPATH="/scripts:$PYTHONPATH"

# bash /scripts/prepdata/all.sh &> /outputs/prepdata.log
# bash /scripts/adjust/all.sh &> /outputs/adjust.log
bash /scripts/metrics/all.sh &> /outputs/metrics.log
# bash /scripts/figures/all.sh &> /outputs/figures.log
