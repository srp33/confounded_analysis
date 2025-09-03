import argparse
import os
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.metrics import make_scorer, roc_auc_score, accuracy_score, mutual_info_score, matthews_corrcoef
from sklearn.model_selection import cross_val_score
from sklearn.feature_selection import mutual_info_classif
from pathlib import Path
from collections import defaultdict
import sys
import time
import random
import numpy as np

# --- Local Imports from utils module ---
from scripts.utils import DataFrameCache, HashCache

from util import repeated_cross_val


def setup_learners():
    # We set random state here, to be updated within repeated_cross_val
    return [
        {"algorithm": RandomForestClassifier, "fit_params": {"n_estimators": 200, "random_state": 0}},
        {"algorithm": HistGradientBoostingClassifier, "fit_params": {"max_iter": 50, "random_state": 0}},
    ]

def mutual_info_shannons(y_true, y_pred):
    """Calculate mutual information in shannons (bits)"""
    return mutual_info_score(y_true, y_pred) / np.log(2)

def setup_metrics():
    return {
        "roc_auc_score": make_scorer(roc_auc_score, response_method=["decision_function", "predict_proba"]),
        "mutual_info_score": make_scorer(mutual_info_shannons),
        "accuracy_score": make_scorer(accuracy_score),
        "matthews_corrcoef": make_scorer(matthews_corrcoef)
    }

def display_metrics(scores_by_metric):
    for metric, values in scores_by_metric.items():
        mean_value = np.mean(values)
        std_value = np.std(values)
        print(f"{metric: <40}  {mean_value:.6f} ± {std_value:.6f}", flush=True)

def initialize_output_file(output_path, debug=False):
    if not os.path.exists(output_path):
        if debug:
            print(f"DEBUG: Initializing output file: {output_path}")
        with open(output_path, 'w') as output_file:
            output_file.write("metric,iteration,classifier,adjuster,dataset,p_column,c_column,c_value,value\n")

def retrieve_and_filter_df(input_path, method, df_cache, index, debug=False):
    df = df_cache.get_dataframe(file_path=input_path)
    if df is None or df.empty:
        return None
    if index is not None:
        # Ensure index is aligned with the dataframe's index
        valid_index = index[index.isin(df.index)]
        if debug:
            print(f"DEBUG: Filtering {len(index)} rows to {len(valid_index)} rows for {input_path}", flush=True)
        df = df.loc[valid_index]
    # If the number of columns is less than 3, return None
    if len(df.columns) < 3:
        return None
    return df

def measure_dataset_method(df_cache, hash_cache, index, learners, method, dataset, prediction_column, conditional_column, c_value, output_path, iterations, input_path, metrics, debug=False):
    df = None # Lazy-load the dataframe only when needed
    for learner in learners:
        classifier_name = learner["algorithm"].__name__.replace("Classifier", "")
        
        key = f"classify|{method}|{dataset}|{prediction_column}|{conditional_column or 'None'}|{c_value}|{classifier_name}"
        description = f"{classifier_name:^25} for {method:^20} on {dataset:>10}, predicting {prediction_column:>25} | {conditional_column}={c_value}"
        
        with hash_cache.check(key, input_path) as should_skip:
            if should_skip:
                print(f"Skipping {description} (no changes)", flush=True)
                continue

            if df is None:
                df = retrieve_and_filter_df(input_path, method, df_cache, index)
                if df is None or df.empty:
                    print(f"Skipping {description} (input not found or empty after filtering)", flush=True)
                    continue

            print(f"Processing {description}", flush=True)
            scores = repeated_cross_val(df, prediction_column, learner, iterations=iterations, n_folds=3, n_jobs=12, metrics=metrics)
            scores_by_metric = {metric: [round(score[metric], 6) for score in scores] for metric in metrics}
            display_metrics(scores_by_metric)

            with open(output_path, 'a') as output_file:
                for metric, values in scores_by_metric.items():
                    for iteration, value in enumerate(values):
                        output_file.write(f"{metric},{iteration},{classifier_name},{method},{dataset},{prediction_column},{conditional_column or 'None'},{c_value},{value}\n")

def process_indexes(df_cache, hash_cache, indexes, input_dir, output_path, prediction_column, conditional_column, write_over, iterations, debug=False):
    dataset = os.path.basename(input_dir)
    
    metrics = setup_metrics()
    learners = setup_learners()

    for index, c_value in indexes:
        # Search by path
        for input_path in input_dir.glob("*.csv"):
            method = input_path.stem
            if method == "unadjusted_t" or method.startswith("."):
                continue
            measure_dataset_method(df_cache, hash_cache, index, learners, method, dataset, prediction_column, conditional_column, c_value, output_path, iterations, input_path, metrics, debug=debug)

def classify_conditionally(df_cache, hash_cache, input_dir, output_path, prediction_column, conditional_column, write_over, iterations, debug=False):
    df_unadjusted_path = input_dir / "unadjusted.csv"
    if not df_unadjusted_path.exists():
        print(f"Warning: unadjusted.csv not found in {input_dir}. Cannot determine conditional values. Skipping.")
        return

    df = df_cache.get_dataframe(file_path=df_unadjusted_path)
    if df.empty:
        print(f"Warning: unadjusted.csv in {input_dir} is empty. Skipping.")
        return
    
    indexes = []
    if conditional_column:
        if conditional_column not in df.columns:
            print(f"Warning: Conditional column '{conditional_column}' not found in {df_unadjusted_path}. Skipping.")
            return
        conditional_values = df[conditional_column].unique()
        for value in conditional_values:
            index = df[df[conditional_column] == value].index
            indexes.append((index, str(value)))
    else:
        indexes.append((None, "None"))

    process_indexes(df_cache, hash_cache, indexes, input_dir, output_path, prediction_column, conditional_column, write_over, iterations, debug=debug)

def main():
    parser = argparse.ArgumentParser(description="Run classification tasks with hash-based caching.")
    parser.add_argument("-i", "--input-dir", type=Path, help="Input directory containing corrected CSVs", required=True)
    parser.add_argument("-o", "--output-path", type=Path, help="Path to output CSV file", required=True)
    parser.add_argument("-p", "--prediction_column", help="Column to predict", required=True)
    parser.add_argument('-c', "--conditional_column", help="Column to condition on", required=False, default=None)
    parser.add_argument("--both_ways", action="store_true", help="Classify both ways (p vs c and c vs p)")
    parser.add_argument("--hash-dir", type=Path, help="Directory to store cache files", required=True)
    parser.add_argument("--write-over", action="store_true", help="Overwrite existing cache files.")
    parser.add_argument("--debug", action="store_true", help="Enable detailed debug prints.")
    
    args = parser.parse_args()
    
    df_cache = DataFrameCache()
    dataset_name = args.input_dir.name
    
    hash_cache = HashCache(
        hash_dir=args.hash_dir,
        # Make sure to store hashes in a specific file, since we want to know if the 
        # input file was changed since *this* script was last run
        cache_filename=f"{dataset_name}_classify.hashes.json",
        write_over=args.write_over,
        debug=True
    )
    
    random.seed(42)

    # Initialize output file once
    initialize_output_file(args.output_path, args.debug)
    
    classify_conditionally(df_cache, hash_cache, args.input_dir, args.output_path, args.prediction_column, args.conditional_column, args.write_over, iterations=10, debug=args.debug)
    
    if args.both_ways and args.conditional_column:
        classify_conditionally(df_cache, hash_cache, args.input_dir, args.output_path, args.conditional_column, args.prediction_column, args.write_over, iterations=10, debug=args.debug)

if __name__ == "__main__":
    main()
