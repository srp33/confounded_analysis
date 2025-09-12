#!/bin/bash

set -e

echo "Running as:"
id

# Set PYTHONPATH for all Python scripts
export PYTHONPATH="/scripts:$PYTHONPATH"

# bash /scripts/prepdata/all.sh &> /outputs/prepdata.log
bash /scripts/adjust/all.sh &> /outputs/adjust.log
# bash /scripts/evaluations/all.sh

