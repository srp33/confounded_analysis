#!/bin/bash

set -e

#bash /scripts/prepdata/all.sh &> /outputs/prepdata.log
# bash /scripts/optimize/all.sh &> /outputs/optimize.log
bash /scripts/adjust/all.sh &> /outputs/adjust.log
bash /scripts/metrics/all.sh &> /outputs/metrics.log
bash /scripts/figures/all.sh &> /outputs/figures.log

# bash /scripts/experiment.sh &> /outputs/experiment.log
# bash /scripts/adjust/profile_adjust.sh &> /outputs/adjust.log