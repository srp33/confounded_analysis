import argparse
import os
import random
from sklearn.svm import SVC
from sklearn.linear_model import SGDClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import RobustScaler
import sys
sys.path.append("/scripts/metrics")
from util import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_path", help="Path to input file", required=True)
parser.add_argument("-o", "--output_path", help="Path to output file", required=True)
parser.add_argument("-b", "--batch_column", help="Batch column", required=True)
parser.add_argument("-p", "--class_column", help="Class column", required=True)
parser.add_argument("-s", "--sample_column", help="Sample identifier column", required=True)
parser.add_argument("-c", "--covariate_columns", help="Covariate columns", required=True)
parser.add_argument("-y", "--variables", help="Variables being evaluated", required=False, default=None)
parser.add_argument("-z", "--variable_values", help="Values of variables being evaluated", required=False, default=None)
args = parser.parse_args()

if args.variables:
    with open(args.output_path, "w") as output_file:
        output_file.write("\t".join(args.variables.split(",")) + "\tdataset\talgorithm\titeration\tcolumn\tvalue\n")

dataset = os.path.basename(args.input_path).replace(".csv", "")
num_cv_iterations = 5
num_folds = 5

LEARNERS = [
    ("Random Forests", RandomForestClassifier(n_estimators = 100)),
    ("Support Vector Machines", SVC(gamma = "auto", kernel = "rbf")),
    ("k-nearest Neighbors", KNeighborsClassifier(n_neighbors = 10))
]

# Setting low_memory=False avoids a warning when you have values (like 4S) that pandas thinks should be numbers.
df = pd.read_csv(args.input_path, low_memory=False)

# Keep only the gene-expression columns
X_cols = df.columns.tolist()
del X_cols[X_cols.index(args.batch_column)]
del X_cols[X_cols.index(args.class_column)]
del X_cols[X_cols.index(args.sample_column)]
for col in args.covariate_columns.split(","):
    del X_cols[X_cols.index(col)]

X = df[X_cols]

num_inf_values = np.isinf(X).values.sum()

if num_inf_values == 0:
    # Scale the X values
    X = RobustScaler().fit(X).transform(X)

    y_batch = df[args.batch_column].to_numpy()
    y_class = df[args.class_column].to_numpy()

with open(args.output_path, "a") as output_file:
    for learner in LEARNERS:
        for i in range(1, num_cv_iterations + 1):
            random.seed(i)

            if num_inf_values == 0:
                print(f"Performing classification for {dataset}, {learner[0]}, iteration {i}.")
                batch_scores = cross_val_score(learner[1], X, y_batch, cv = num_folds, scoring = "roc_auc", n_jobs = 1)
                class_scores = cross_val_score(learner[1], X, y_class, cv = num_folds, scoring = "roc_auc", n_jobs = 1)

                batch_score = float(np.mean(batch_scores))
                class_score = float(np.mean(class_scores))
            else:
                print(f"NOT performing classification for {dataset}, {learner[0]}, iteration {i} because there are inf values.")
                batch_score = "NA"
                class_score = "NA"

            if args.variable_values:
                variable_values = args.variable_values.split(",")
            else:
                variable_values = ["NA" for x in args.variables.split(",")]

            output_file.write("\t".join(variable_values) + f"\t{dataset}\t{learner[0]}\t{i}\tBatch\t{batch_score}\n")
            output_file.write("\t".join(variable_values) + f"\t{dataset}\t{learner[0]}\t{i}\tClass\t{class_score}\n")
