import argparse
import os
from sklearn.svm import SVC
from sklearn.linear_model import SGDClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
import sys
sys.path.append("/scripts/metrics")
from util import *

#import os.path
#from os import path

#import time
#import random

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input-path", help="Path to input file", required=True)
parser.add_argument("-o", "--output-path", help="Path to output file", required=True)
parser.add_argument("-b", "--batch-column", help="Batch column", required=True)
parser.add_argument("-p", "--pred-column", help="Prediction column", required=True)
parser.add_argument("-y", "--variables", help="Variables being evaluated", required=True)
parser.add_argument("-z", "--variable_values", help="Values of variables being evaluated", required=True)
args = parser.parse_args()

cache = DataFrameCache()

print(args.output_path)
if not os.path.exists(args.output_path):
    with open(args.output_path, "w") as output_file:
        output_file.write("\t".join(args.variables.split(",")) + "\tdataset\talgorithm\tmetric\tvalue\n")

#results = []
#random.seed()
#
#dataset = os.path.basename(args.input_dir)
#
#TODO: Move random_state later
#LEARNERS = [
#    (RandomForestClassifier, {"n_estimators": 100, "random_state": 0}),
#    (SVC, {"gamma": "auto", "random_state": 0, "kernel": "rbf"}),
#    (KNeighborsClassifier, {})
#]

#for method in ["unadjusted", "scaled", "combat", "confounded"]:
#    df = cache.get_dataframe(args.input_dir + "/" + method + ".csv")
#
#    for learner in LEARNERS:
#        classifier_name = str(learner[0]).split("'")[1].split(".")[-1].replace("Classifier", "")
#
#        print("Performing classification for {}, {}, {}, and {}.".format(dataset, method, args.column, classifier_name), flush=True)
#        for score in cross_validate(df, args.column, learner, iterations=10, folds=3, n_jobs=12):
#           results.append([classifier_name, method, dataset, str(score)])
#
#with open(args.output_path, 'a') as output_file:
#    for line in results:
#        output_file.write(",".join(line) + "\n")
