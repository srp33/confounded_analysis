#!/bin/bash

# Set PYTHONPATH for all Python scripts
export PYTHONPATH="/scripts:$PYTHONPATH"

bash /scripts/evaluations/classify_er_mixed_datasets/hist_gradient_er.sh &> /outputs/hist_gradient_er.log