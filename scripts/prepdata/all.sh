#!/bin/bash

set -e

bash /scripts/prepdata/simulate_expression.sh
bash /scripts/prepdata/bladderbatch.sh
bash /scripts/prepdata/gse37199.sh
bash /scripts/prepdata/tcga.sh
