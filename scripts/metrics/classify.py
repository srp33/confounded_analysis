import argparse
import os
import pandas as pd
from sklearn.svm import SVC
from sklearn.linear_model import SGDClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import robust_scale
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
import os.path
from os import path

import sys
import time
import random
from util import DataFrameCache, get_dataset_path_dict, no_extension

def cross_validate(df, predict_column, learner):
    meta_cols = list(df.select_dtypes(include=['object', 'int']).columns)
    X = robust_scale(df.drop(meta_cols, axis="columns"))
    y = df[predict_column]

    #classes = []
    #for element in y:
    #    if (element not in classes) : 
    #        classes.append(element)
    #y = label_binarize(y, classes = classes)
    #print(y)

    scoring_metric = "accuracy"
    n_jobs = 12

    scores = []
    for i in range(iterations):
        fit_params = learner[1]
        if "random_state" in fit_params:
            fit_params["random_state"] = i

        estimator = learner[0](**fit_params)

        kfold = StratifiedKFold(n_splits=folds, shuffle=True, random_state=random.randint(1,100))
        iter_scores = list(cross_val_score(estimator, X, y, scoring=scoring_metric, cv=kfold, n_jobs=n_jobs))

        scores.append(sum(iter_scores) / len(iter_scores))

    return scores

def baseline(df, column):
    return df[column].value_counts().max() / len(df)

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input-dir", help="Input directory", required=True)
parser.add_argument("-o", "--output-path", help="Path to output file", required=True)
parser.add_argument("-c", "--column", help="Prediction column", required=True)
args = parser.parse_args()

iterations = 10
# It makes sense to use few folds because there are few samples for some of the classes
folds = 3

cache = DataFrameCache()

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

#results.append(["baseline", "NA", dataset, str(baseline(unadj, args.column))])

for method in ["unadjusted", "scaled", "combat", "confounded"]:
    df = cache.get_dataframe(args.input_dir + "/" + method + ".csv")
    for learner in LEARNERS:
        classifier_name = str(learner[0]).split("'")[1].split(".")[-1].replace("Classifier", "")

        print("Performing classification for {}, {}, {}, and {}.".format(dataset, method, args.column, classifier_name), flush=True)
        for score in cross_validate(df, args.column, learner):
           results.append([classifier_name, method, dataset, str(score)])

with open(args.output_path, 'a') as output_file:
    for line in results:
        output_file.write(",".join(line) + "\n")
