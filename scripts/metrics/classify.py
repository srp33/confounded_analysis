import argparse
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NeighborhoodComponentsAnalysis
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import make_scorer
from sklearn.metrics import roc_auc_score
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
    {"algorithm": RandomForestClassifier, "fit_params": {"n_estimators": 200, "random_state": 0}, 
     "transform": NeighborhoodComponentsAnalysis, "transform_params": {"n_components": 50}},
    {"algorithm": GradientBoostingClassifier, "fit_params": {"n_estimators": 50, "random_state": 0}},
    {"algorithm": HistGradientBoostingClassifier, "fit_params": {"max_iter": 50, "random_state": 0}},
    {"algorithm": KNeighborsClassifier, "fit_params": {"n_neighbors": 5, "metric": "cosine"}}
]

if not os.path.exists(args.output_path):
    with open(args.output_path, "w") as output_file:
        output_file.write("metric,adjuster,dataset,column,value\n")

results = []
random.seed()

dataset = os.path.basename(args.input_dir)


# For mutual_info_score, we want to return the score in shannons (bits)
def mutual_info_shannons(y_true, y_pred):
    """Calculate mutual information in shannons (bits)"""
    return mutual_info_score(y_true, y_pred) / np.log(2)



metrics = [
    make_scorer(roc_auc_score, needs_proba=True, multi_class='ovr', average='macro'),
    make_scorer(mutual_info_shannons, needs_proba=False)
]



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
            for metric in metrics:
                score_for_metric = score[metric]

                results.append([classifier_name, method, dataset, args.column, metric, str(score_for_metric)])

with open(args.output_path, 'a') as output_file:
    for line in results:
        output_file.write(",".join(line) + "\n")
