import argparse
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier, NeighborhoodComponentsAnalysis
from sklearn.ensemble import GradientBoostingClassifier, HistGradientBoostingClassifier
from sklearn.metrics import make_scorer, roc_auc_score, accuracy_score, log_loss, mutual_info_score
from sklearn.model_selection import cross_val_score


import os.path
from os import path

import sys
import time
import random
from util import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input-dir", help="Input directory", required=True)
parser.add_argument("-o", "--output-path", help="Path to output file", required=True)
parser.add_argument("-c", "--column", help="Prediction column", required=True)
parser.add_argument('--class0', nargs='+', required=False, help='Values to map to class 0')
parser.add_argument('--class1', nargs='+', required=False, help='Values to map to class 1')
args = parser.parse_args()

cache = DataFrameCache()

nca = NeighborhoodComponentsAnalysis(n_components=100, random_state=42)

LEARNERS = [  # Random state is updated within repeated_cross_val
    {"algorithm": RandomForestClassifier, "fit_params": {"n_estimators": 200, "random_state": 0}},
    {"algorithm": HistGradientBoostingClassifier, "fit_params": {"max_iter": 50, "random_state": 0}},
]

if not os.path.exists(args.output_path):
    with open(args.output_path, "w") as output_file:
        output_file.write("metric,adjuster,dataset,column,score,value\n")

results = []
random.seed()

dataset = os.path.basename(args.input_dir)


# For mutual_info_score, we want to return the score in shannons (bits)
def mutual_info_shannons(y_true, y_pred):
    """Calculate mutual information in shannons (bits)"""
    return mutual_info_score(y_true, y_pred) / np.log(2)


def mutual_info_proba_shannons(y_true, y_pred_proba):
    """
    Calculate mutual information between true classes (categorical)
    and predicted probabilities (continuous) in shannons (bits).
    """
    # For binary classification, y_pred_proba is usually a 2D array: [[prob_class_0, prob_class_1], ...]
    # We can use either column as our continuous variable.
    continuous_output = y_pred_proba[:, 0]

    # mutual_info_classif expects X to be 2D, even for a single feature.
    # So, reshape continuous_output.
    mi_value = mutual_info_classif(continuous_output.reshape(-1, 1), y_true, random_state=42)[0]

    # Convert nats (default for mutual_info_classif) to shannons (bits)
    # 1 nat = 1 / log(2) bits
    return mi_value / np.log(2)


metrics = {
    "roc_auc_score": make_scorer(roc_auc_score, response_method=["decision_function", "predict_proba"]),
    "mutual_info_score": make_scorer(mutual_info_shannons),
    "accuracy_score": make_scorer(accuracy_score),
    "log_loss": make_scorer(log_loss, response_method=["predict_proba"], greater_is_better=False)
}



for method in ["combat", "combat_target", "limma_target", "unadjusted", "min_mean", "tampor", "limma"]:
    df = cache.get_dataframe(args.input_dir + "/" + method + ".csv")

    if args.class0 and args.class1:
        # Map classes and flag invalid entries
        valid_classes = set(args.class0 + args.class1)
        df[args.column] = df[args.column].astype(str)

        invalid = df[~df[args.column].isin(valid_classes)]
        if not invalid.empty:
            print(f"Unmapped classes found: {invalid[args.column].unique()}")
        
        df[args.column] = df[args.column].apply(
            lambda x: 0 if x in args.class0 else 1 if x in args.class1 else np.nan
        ).dropna()  # Remove rows with unmapped classes
    
    for learner in LEARNERS:
        classifier_name = str(learner["algorithm"]).split("'")[1].split(".")[-1].replace("Classifier", "")

        print("Performing classification for {}, {}, {}, and {}.".format(dataset, method, args.column, classifier_name), flush=True)
        for score in repeated_cross_val(df, args.column, learner, iterations=10, n_folds=3, n_jobs=12, metrics=metrics):
            print(f"Score: {score}", flush=True)
            for metric in metrics:
                score_for_metric = score[metric]

                results.append([classifier_name, method, dataset, args.column, metric, str(score_for_metric)])

with open(args.output_path, 'a') as output_file:
    for line in results:
        output_file.write(",".join(line) + "\n")
