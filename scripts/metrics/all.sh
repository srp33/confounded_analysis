#!/bin/bash

set -e

# python /scripts/metrics/grab_esr1.py

# bash /scripts/metrics/feature_importance.sh

# python /scripts/metrics/esr1_analysis.py


# bash /scripts/metrics/mutual_info.sh
# bash /scripts/metrics/classify.sh
# bash /scripts/metrics/classify_combined.sh
bash /scripts/metrics/hist_gradient_er.sh
# bash /scripts/metrics/reduce.sh
# bash /scripts/metrics/mse.sh
# bash /scripts/metrics/mmd.sh

# python /scripts/metrics/eval_refine.py