#!/bin/bash

set -e

# python /scripts/metrics/hist_gradient_er_classification.py

# bash /scripts/metrics/mutual_info.sh
bash /scripts/metrics/reduce.sh
bash /scripts/metrics/classify.sh
# bash /scripts/metrics/mse.sh
# bash /scripts/metrics/mmd.sh

# python /scripts/metrics/eval_refine.py