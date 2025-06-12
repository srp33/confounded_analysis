#!/bin/bash

set -e

#bash /scripts/metrics/mse.sh
#bash /scripts/metrics/mmd.sh
bash /scripts/metrics/classify.sh
# Mutual info depends on classify
bash /scripts/metrics/mutual_info.sh
