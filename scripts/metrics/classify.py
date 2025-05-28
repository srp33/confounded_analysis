import argparse
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier

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

LEARNERS = [
    (RandomForestClassifier, {"n_estimators": 100, "random_state": 0}), # Random state is updated within cross_validate
    (MLPClassifier, {"alpha": 0.01, "max_iter": 1000, "random_state": 0}),
    (KNeighborsClassifier, {"n_neighbors": 7, "weights": "distance", "metric": "cosine"}),
]

if not os.path.exists(args.output_path):
    with open(args.output_path, "w") as output_file:
        output_file.write("metric,adjuster,dataset,column,value\n")

results = []
random.seed()

dataset = os.path.basename(args.input_dir)

for method in ["unadjusted", "min_mean", "combat", "tampor"]:
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
        classifier_name = str(learner[0]).split("'")[1].split(".")[-1].replace("Classifier", "")

        print("Performing classification for {}, {}, {}, and {}.".format(dataset, method, args.column, classifier_name), flush=True)
        for score in cross_validate(df, args.column, learner, iterations=10, folds=3, n_jobs=12):
           results.append([classifier_name, method, dataset, args.column, str(score)])

with open(args.output_path, 'a') as output_file:
    for line in results:
        output_file.write(",".join(line) + "\n")
