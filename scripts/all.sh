#!/bin/bash

set -e

#bash /scripts/prepdata/all.sh &> /outputs/prepdata.log
#bash /scripts/optimize/all.sh &> /outputs/optimize.log
bash /scripts/optimize/all.sh
#TODO: Add https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24080
#bash /scripts/adjust/all.sh &> /outputs/adjust.log
#bash /scripts/metrics/all.sh &> /outputs/metrics.log
#bash /scripts/figures/all.sh &> /outputs/figures.log
