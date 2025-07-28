import argparse
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier, NeighborhoodComponentsAnalysis
from sklearn.ensemble import GradientBoostingClassifier, HistGradientBoostingClassifier
from sklearn.metrics import make_scorer, roc_auc_score, accuracy_score, log_loss, mutual_info_score
from sklearn.model_selection import cross_val_score
from sklearn.feature_selection import mutual_info_classif
from pathlib import Path
from collections import defaultdict


import os.path
from os import path

import sys
import time
import random
from util import *


def setup_learners():
    return [  # Random state is updated within repeated_cross_val
    {"algorithm": RandomForestClassifier, "fit_params": {"n_estimators": 200, "random_state": 0}},
    {"algorithm": HistGradientBoostingClassifier, "fit_params": {"max_iter": 50, "random_state": 0}},
]


def setup_metrics():
    return {
    "roc_auc_score": make_scorer(roc_auc_score, response_method=["decision_function", "predict_proba"]),
    "mutual_info_score": make_scorer(mutual_info_shannons),
    # "mutual_info_proba_classif": make_scorer(mutual_info_proba_classif, response_method=["predict_proba"]),
    # "mutual_info_proba_sum": make_scorer(mutual_info_proba_sum, response_method=["predict_proba"]),
    # "determinant_based_mutual_information": make_scorer(determinant_based_mutual_information, response_method=["predict_proba"]),
    "accuracy_score": make_scorer(accuracy_score),
    # "log_loss": make_scorer(one_minus_log_loss, response_method=["predict_proba"])
}


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


def mine_previous_results(file_path, iterations):
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

    previous_results = defaultdict(dict)
    for line in [l.strip() for l in lines if l.strip()]:
        metric, iteration, classifier, adjuster, dataset, p_column, c_column, c_value, value = line.split(",")
        key = (classifier, adjuster, dataset, p_column, c_column, c_value)
        if key not in previous_results:
            previous_results[key] = {}
        if metric not in previous_results[key]:
            previous_results[key][metric] = [-1] * iterations
        previous_results[key][metric][int(iteration)] = value

    return previous_results


def display_metrics(scores_by_metric):
    for metric, values in scores_by_metric.items():
        mean_value = np.mean(values)
        std_value = np.std(values)
        print(f"{metric: <40}  {mean_value:.6f} ± {std_value:.6f}", flush=True)


def initialize_output_file(output_path, write_over):
    if write_over or not path.exists(output_path):
        with open(output_path, 'w') as output_file:
            output_file.write("metric,iteration,classifier,adjuster,dataset,p_column,c_column,c_value,value\n")



def retrieve_and_filter_df(input_path, method, cache, index, classifier_name, dataset):
    df = cache.get_dataframe(input_path)
    print(f"Getting {input_path}", flush=True)
    print(f"Total rows: {len(df)}", flush=True)
    if df is None or len(df) == 0:
        return None
    if index is not None:
        df = df.loc[index]
    return df



def measure_dataset_method(cache, index, learners, method, dataset, prediction_column, conditional_column, c_value, results, output_path, write_over, iterations, input_path, metrics):
    df = None
    for learner in learners:
        classifier_name = str(learner["algorithm"]).split("'")[1].split(".")[-1].replace("Classifier", "")
        key = (classifier_name, method, dataset, prediction_column, str(conditional_column), c_value)
        description = f"{classifier_name:^25} for {method:^20} on {dataset:>10}, predicting {prediction_column:>15} | {conditional_column}={c_value}"

        if key in results and not write_over:
            print(f"Skipping {description}", flush=True)
            continue

        if df is None:
            df = retrieve_and_filter_df(input_path, method, cache, index, classifier_name, dataset)
            if df is None:
                print(f"Skipping {description}", flush=True)
                continue

        print(f"Performing classification using {classifier_name} for {method} on {dataset}, predicting {prediction_column} | {conditional_column}={c_value}", flush=True)
        scores = repeated_cross_val(df, prediction_column, learner, iterations=iterations, n_folds=3, n_jobs=12, metrics=metrics)

        # Reduce to 6 decimal places
        scores_by_metric = {metric: [round(score[metric], 6) for score in scores] for metric in metrics}

        # Print mean and std of each metric
        display_metrics(scores_by_metric)

        # Transpose structure, write to file
        with open(output_path, 'a') as output_file:
            for metric, values in scores_by_metric.items():
                results[key][metric] = values
                for iteration, value in enumerate(values):
                    output_file.write(f"{metric},{iteration},{classifier_name},{method},{dataset},{prediction_column},{conditional_column},{c_value},{value}\n")



def process_indexes(cache, indexes, input_dir, output_path, prediction_column, conditional_column, write_over, iterations):
    dataset = os.path.basename(input_dir)
    results = mine_previous_results(output_path, iterations)
    print(f"{len(results)} previous results found in {output_path}", flush=True)

    metrics = setup_metrics()
    learners = setup_learners()
    methods = ["min_mean", "combat", "combat_target", "npn", "limma", "limma_target", "unadjusted", "tampor", "quantile", "autoclass", "icvae", "seurat_scaling", "seurat_integration", "fastMNN", "liger", "wasserstein", "monotonic", "non_monotonic"]

    for index, c_value in indexes:
        for method in methods:
            input_path = input_dir / f"{method}.csv"
            measure_dataset_method(cache, index, learners, method, dataset, prediction_column, conditional_column, c_value, results, output_path, write_over, iterations, input_path, metrics)



def classify_conditionally(cache, input_dir, output_path, prediction_column, conditional_column, write_over, iterations):
    df = cache.get_dataframe(input_dir / "unadjusted.csv")
    print(f"Length of unadjusted: {len(df)}", flush=True)

    indexes = []
    if conditional_column:
        conditional_values = df[conditional_column].unique()
        for value in conditional_values:
            index = df[conditional_column] == value
            indexes.append((index, str(value)))
            print(f"Length of {conditional_column}={value}: {len(df[df[conditional_column] == value])}", flush=True)
    else:
        # None represents "all rows"
        indexes.append((None, "None"))

    process_indexes(cache, indexes, input_dir, output_path, prediction_column, conditional_column, write_over, iterations)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-dir", type=Path, help="Input directory", required=True)
    parser.add_argument("-o", "--output-path", help="Path to output file", required=True)
    parser.add_argument("-p", "--prediction_column", help="Prediction column", required=True)
    parser.add_argument('-c', "--conditional_column", help="Conditional column", required=False)
    parser.add_argument("--both_ways", action="store_true", help="Classify both ways", required=False)
    parser.add_argument("-w", "--write-over", action="store_true", help="Redo computations even if done before, writing over the previous results")
    args = parser.parse_args()
    cache = DataFrameCache()
    random.seed()

    initialize_output_file(args.output_path, args.write_over)
    classify_conditionally(cache, args.input_dir, args.output_path, args.prediction_column, args.conditional_column, args.write_over, iterations=10)
    if args.both_ways:
        classify_conditionally(cache, args.input_dir, args.output_path, args.conditional_column, args.prediction_column, args.write_over, iterations=10)
    

if __name__ == "__main__":
    main()
