import argparse
import os
from sklearn.svm import SVC
from sklearn.linear_model import SGDClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression

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
args = parser.parse_args()

cache = DataFrameCache()

#TODO: Move random_state later
LEARNERS = [
    (RandomForestClassifier, {"n_estimators": 100, "random_state": 0}),
    (SVC, {"gamma": "auto", "random_state": 0, "kernel": "rbf"}),
    (KNeighborsClassifier, {})
]

if not os.path.exists(args.output_path):
    with open(args.output_path, "w") as output_file:
        output_file.write("metric,adjuster,dataset,value\n")

results = []
random.seed()

dataset = os.path.basename(args.input_dir)

for method in ["unadjusted", "scaled", "combat", "confounded"]:
    df = cache.get_dataframe(args.input_dir + "/" + method + ".csv")

    for learner in LEARNERS:
        classifier_name = str(learner[0]).split("'")[1].split(".")[-1].replace("Classifier", "")

        print("Performing classification for {}, {}, {}, and {}.".format(dataset, method, args.column, classifier_name), flush=True)
        for score in cross_validate(df, args.column, learner, iterations=10, folds=3, n_jobs=12):
           results.append([classifier_name, method, dataset, str(score)])

with open(args.output_path, 'a') as output_file:
    for line in results:
        output_file.write(",".join(line) + "\n")
