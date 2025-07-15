import argparse
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier, NeighborhoodComponentsAnalysis
from sklearn.ensemble import GradientBoostingClassifier, HistGradientBoostingClassifier
from sklearn.metrics import make_scorer, roc_auc_score, accuracy_score, log_loss, mutual_info_score
from sklearn.model_selection import cross_val_score
from sklearn.feature_selection import mutual_info_classif


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
parser.add_argument("-w", "--write-over", action="store_true", help="Redo computations even if done before, writing over the previous results")
args = parser.parse_args()

cache = DataFrameCache()

nca = NeighborhoodComponentsAnalysis(n_components=100, random_state=42)

LEARNERS = [  # Random state is updated within repeated_cross_val
    {"algorithm": RandomForestClassifier, "fit_params": {"n_estimators": 200, "random_state": 0}},
    {"algorithm": HistGradientBoostingClassifier, "fit_params": {"max_iter": 50, "random_state": 0}},
]

ITERATIONS = 10

results = []
random.seed()

dataset = os.path.basename(args.input_dir)


# For mutual_info_score, we want to return the score in shannons (bits)
def mutual_info_shannons(y_true, y_pred):
    """Calculate mutual information in shannons (bits)"""
    return mutual_info_score(y_true, y_pred) / np.log(2)


def mutual_info_proba_classif(y_true, y_pred_proba):
    """
    Calculate mutual information between true classes (categorical)
    and predicted probabilities (continuous) in shannons (bits), using mutual_info_classif.
    """
    # For binary classification, y_pred_proba is a 1D array: [n_samples]
    # For multiclass classification, y_pred_proba is a 2D array: [n_samples, n_classes]
    continuous_output = y_pred_proba

    # mutual_info_classif expects X to be 2D, even for a single feature.
    # So, reshape continuous_output.
    mi_value = mutual_info_classif(continuous_output.reshape(-1, 1), y_true, random_state=42)[0]

    # Convert nats (default for mutual_info_classif) to shannons (bits)
    # 1 nat = 1 / log(2) bits
    return mi_value / np.log(2)


def mutual_info_proba_sum(y_true, y_pred_proba, classes=[]):
    """
    Calculate mutual information by summing the probabilities across samples, for each true class.
    This results in a contingency table (2x2 for binary classification), which is then used to compute mutual information.
    This is similar to mutual_info_score but sums p(y_pred_proba | y_true) across samples, instead of 1-0 encoding the probabilities.
    """
    # For binary classification, y_pred_proba is a 1D array: [n_samples]
    # For multiclass classification, y_pred_proba is a 2D array: [n_samples, n_classes]
    binary = len(y_pred_proba.shape) == 1

    if not classes:
        classes = sorted(np.unique(y_true).tolist())

    # Create a contingency table
    contingency_table = {
        cls: np.zeros(len(classes)) for cls in classes
    }

    if binary:
        # Add 1-pred as a column to y_pred_proba
        y_pred_proba = np.column_stack((y_pred_proba, 1 - y_pred_proba))


    for true_class, proba in zip(y_true, y_pred_proba):
        contingency_table[true_class] += proba

    # Convert table to np array
    # Contingency table accepts ints only. Multiplying by a large number prevents rounding errors.
    contingency_table = np.array(list(contingency_table.values())) * 10**6

    # Calculate mutual information from the contingency table
    return mutual_info_score(None, None, contingency=contingency_table) / np.log(2)


def one_minus_log_loss(y_true, y_pred_proba):
    """
    Calculate one minus log loss.
    This is useful for maximizing the score, as log loss is minimized.
    """
    return 1 - log_loss(y_true, y_pred_proba)


def determinant_based_mutual_information(y_true, y_pred_proba):
    """
    Calculate mutual information using det(P), where P is the joint probability distribution of y_true and y_pred_proba.
    """
    classes = sorted(np.unique(y_true).tolist())

    # Convert y_true to one-hot encoding
    one_hot = np.zeros((len(y_true), len(classes)))
    for i, cls in enumerate(classes):
        one_hot[:, i] = (y_true == cls).astype(int)
    y_true = one_hot

    y_pred_proba = np.asarray(y_pred_proba)
    
    binary = len(y_pred_proba.shape) == 1
    if binary:
        # For binary classification, y_pred_proba is a 1D array: [n_samples]
        # For multiclass classification, y_pred_proba is a 2D array: [n_samples, n_classes]
        y_pred_proba = np.column_stack((y_pred_proba, 1 - y_pred_proba))

    contingency_table = y_pred_proba.T @ y_true
    rel_freq_table = contingency_table / contingency_table.sum()

    # Calculate the determinant of the contingency table
    det = np.linalg.det(rel_freq_table + 1e-10)  # Add a small value to avoid singular matrix issues
    return abs(det) * 4 # Max is 0.25 (max of x*1-x) because the probabilities must sum to 1. [0.5, 0] x [0, 0.5] = 0.25


def mine_previous_results(file_path):
    """
    Read the output file and return a dictionary of previous results.
    The keys are tuples of (metric, iteration, adjuster, dataset, column).
    The values are the scores.
    """
    if not path.exists(file_path):
        return {}
    with open(file_path, 'r') as output_file:
        lines = output_file.readlines()

    # Skip the header line
    lines = lines[1:]

    previous_results = {}
    for line in [l.strip() for l in lines if l.strip()]:
        metric, iteration, classifier, adjuster, dataset, column, value = line.split(",")
        key = (classifier, adjuster, dataset, column)
        if key not in previous_results:
            previous_results[key] = {}
        if metric not in previous_results[key]:
            previous_results[key][metric] = [-1] * ITERATIONS
        previous_results[key][metric][int(iteration)] = value

    return previous_results


metrics = {
    "roc_auc_score": make_scorer(roc_auc_score, response_method=["decision_function", "predict_proba"]),
    "mutual_info_score": make_scorer(mutual_info_shannons),
    "mutual_info_proba_classif": make_scorer(mutual_info_proba_classif, response_method=["predict_proba"]),
    "mutual_info_proba_sum": make_scorer(mutual_info_proba_sum, response_method=["predict_proba"]),
    "determinant_based_mutual_information": make_scorer(determinant_based_mutual_information, response_method=["predict_proba"]),
    "accuracy_score": make_scorer(accuracy_score),
    "log_loss": make_scorer(one_minus_log_loss, response_method=["predict_proba"])
}

def binarize_column(df, column, class0, class1):
    """
    Binarize the specified column in the DataFrame.
    Map class0 to 0, class1 to 1, and drop rows with unmapped classes.
    """
    valid_classes = set(class0 + class1)
    df[column] = df[column].astype(str)

    # Flag invalid entries
    invalid = df[~df[column].isin(valid_classes)]
    if not invalid.empty:
        print(f"Unmapped classes found: {invalid[column].unique()}")

    # Map classes
    df[column] = df[column].apply(
        lambda x: 0 if x in class0 else 1 if x in class1 else np.nan
    ).dropna()  # Remove rows with unmapped classes

    return df


def display_metrics(scores_by_metric):
    for metric, values in scores_by_metric.items():
        mean_value = np.mean(values)
        std_value = np.std(values)
        print(f"{metric: <40}  {mean_value:.6f} ± {std_value:.6f}", flush=True)

# Read output file to check if we need to redo computations
if args.write_over or not path.exists(args.output_path):
    with open(args.output_path, 'w') as output_file:
        output_file.write("metric,iteration,classifier,adjuster,dataset,column,value\n")

results = mine_previous_results(args.output_path)
print(f"{len(results)} previous results found in {args.output_path}", flush=True)

for method in ["combat", "combat_target", "limma", "limma_target", "unadjusted", "tampor", "quantile", "autoclass", "icvae", "seurat_scaling", "seurat_integration", "fastMNN", "liger", "wasserstein"]:
    df = None
    
    for learner in LEARNERS:
        classifier_name = str(learner["algorithm"]).split("'")[1].split(".")[-1].replace("Classifier", "")
        key = (classifier_name, method, dataset, args.column)
        if not args.write_over and key in results:
            print(f"Skipping {classifier_name} for {method}, {dataset}, {args.column} as results already exist.", flush=True)
            continue
        if key not in results:
            results[key] = {}

        if df is None:
            df = cache.get_dataframe(args.input_dir + "/" + method + ".csv")
            if args.class0 and args.class1:
                df = binarize_column(df, args.column, args.class0, args.class1)

        print(f"Performing classification using {classifier_name} for {method}, {dataset}, {args.column}", flush=True)
        scores = repeated_cross_val(df, args.column, learner, iterations=ITERATIONS, n_folds=3, n_jobs=12, metrics=metrics)

        # Reduce to 6 decimal places, and transpose
        scores_by_metric = {metric: [round(score[metric], 6) for score in scores] for metric in metrics}
        for metric, values in scores_by_metric.items():
            results[key][metric] = values

        # Print mean and std of each metric
        display_metrics(scores_by_metric)

        # Write results to output file
        with open(args.output_path, 'a') as output_file:
            for metric, values in scores_by_metric.items():
                for iteration, value in enumerate(values):
                    output_file.write(f"{metric},{iteration},{classifier_name},{method},{dataset},{args.column},{value}\n")
        

