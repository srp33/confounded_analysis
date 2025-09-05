#!/bin/bash

set -e

# python /scripts/evaluations/esr1/esr1_analysis.py &> /outputs/esr1_analysis.log

# python /scripts/evaluations/quick_classify_refine_datasets/eval_refine.py &> /outputs/eval_refine.log

# bash /scripts/evaluations/lassifier_feature_importance/feature_importance.sh &> /outputs/feature_importance.log

# bash /scripts/evaluations/reduce_data_for_viewing/reduce.sh &> /outputs/reduce.log
# bash scripts/evaluations/reduce_data_for_viewing/plot_reduced.sh &> /outputs/plot_reduce.log

# bash /scripts/evaluations/small_evals/mutual_info.sh &> /outputs/mutual_info.log
# bash /scripts/evaluations/small_evals/mse.sh &> /outputs/mse.log
# bash /scripts/evaluations/small_evals/mmd.sh &> /outputs/mmd.log
# bash scripts/evaluations/small_evals/mse_mmd_classification.sh &> /outputs/mse_mmd_classification.log

# bash /scripts/evaluations/classify_batch_bio_within_dataset/classify.sh &> /outputs/classify.log
# bash /scripts/evaluations/classify_batch_bio_within_dataset/classify_combined.sh &> /outputs/classify_combined.log
# bash scripts/evaluations/classify_batch_bio_within_dataset/classification_figures.sh &> /outputs/classification_figures.log

bash /scripts/evaluations/classify_er_mixed_datasets/hist_gradient_er.sh &> /outputs/hist_gradient_er.log
Rscript scripts/evaluations/classify_er_mixed_datasets/er_classification_plots_single.R &> /outputs/er_classification_plots_single.log






