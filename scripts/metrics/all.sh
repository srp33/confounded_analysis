#!/bin/bash

set -e

# bash /scripts/metrics/mutual_info.sh
# bash /scripts/metrics/mse.sh
# bash /scripts/metrics/mmd.sh
bash /scripts/metrics/classify.sh

# python /scripts/metrics/eval_refine.py