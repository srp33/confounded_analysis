#!/bin/bash

set -e

# bash /scripts/adjust/autoclass.sh
# bash /scripts/adjust/icvae.sh
# bash /scripts/adjust/vfae.sh
# bash /scripts/adjust/wasserstein.sh

# bash /scripts/adjust/adjustR_data.sh
# bash /scripts/adjust/adjustR_individual_prep.sh

echo "🔗 Generating all dataset combinations with caching (gmm and npn files)..."
# python3 scripts/prepdata/generate_all_combinations.py \
#     --data-dir /data/gold \
#     --csv-files \
#         "npn_global.csv" \
#         "gmm_global.csv" \
#         "gmm_npn_global.csv" \
#         "gmm_npn_unit_std_global.csv" \
#         "gmm_scale_separate_global.csv" \
#     --debug

bash /scripts/adjust/adjustR_combined_data.sh


