import argparse
import os
import random
from sklearn.svm import SVC
from sklearn.linear_model import SGDClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
import sys
sys.path.append("/scripts/metrics")
from util import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_path", help="Path to input file", required=True)
parser.add_argument("-o", "--output_path", help="Path to output file", required=True)
parser.add_argument("-b", "--batch_column", help="Batch column", required=True)
parser.add_argument("-p", "--class_column", help="Class column", required=True)
parser.add_argument("-y", "--variables", help="Variables being evaluated", required=False, default=None)
parser.add_argument("-z", "--variable_values", help="Values of variables being evaluated", required=False, default=None)
parser.add_argument("-s", "--scale_numerics", help="Whether to scale the numeric data", required=False, default=False)
args = parser.parse_args()

if args.variables:
    with open(args.output_path, "w") as output_file:
        output_file.write("\t".join(args.variables.split(",")) + "\tdataset\talgorithm\tcolumn\tvalue\n")

random.seed(0)
cache = DataFrameCache()
dataset = os.path.basename(args.input_path).replace(".csv", "")
num_cv_iterations = 2
num_folds = 10
num_jobs = 1

LEARNERS = [
    (RandomForestClassifier, {"n_estimators": 100, "random_state": None}),
    (SVC, {"gamma": "auto", "kernel": "rbf", "random_state": None}),
    (KNeighborsClassifier, {"n_neighbors": 10})
]

df = cache.get_dataframe(args.input_path)

with open(args.output_path, "a") as output_file:
    for learner in LEARNERS:
        classifier_name = learner[0].__name__

        print(f"Performing classification for {dataset}, {classifier_name}.")
        batch_scores = cross_validate(df, args.batch_column, learner, iterations=num_cv_iterations, folds=num_folds, n_jobs=num_jobs, scale_numerics=args.scale_numerics == "True")
        class_scores = cross_validate(df, args.class_column, learner, iterations=num_cv_iterations, folds=num_folds, n_jobs=num_jobs, scale_numerics=args.scale_numerics == "True")

        for i in range(num_cv_iterations):
            if args.variable_values:
                variable_values = args.variable_values.split(",")
            else:
                variable_values = ["NA" for x in args.variables.split(",")]

            output_file.write("\t".join(variable_values) + f"\t{dataset}\t{classifier_name}\tBatch\t{batch_scores[i]}\n")
            output_file.write("\t".join(variable_values) + f"\t{dataset}\t{classifier_name}\tClass\t{class_scores[i]}\n")
