import argparse
import os
import sys
import traceback
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import roc_auc_score, confusion_matrix
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.base import clone
from joblib import Parallel, delayed
from filelock import FileLock

# Add the parent directory (scripts) to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def safe_write_to_csv(df, output_file):
    lock_path = output_file + '.lock'
    with FileLock(lock_path):
        df.to_csv(output_file, mode='a', header=False, index=False, float_format='%.4f')

def print_now(*args, **kwargs):
    """Print a message to the console with flushing."""
    print(*args, flush=True, **kwargs)

def calculate_metrics(y_true, y_pred, y_proba):
    """Calculate a standard set of classification metrics."""
    # Return NaN for metrics that fail if only one class is present.
    if len(pd.unique(y_true)) < 2:
        return {
            'ROC AUC': np.nan,
            'True Negative': np.nan,
            'False Positive': np.nan,
            'False Negative': np.nan,
            'True Positive': np.nan
        }

    # Use confusion matrix for robust calculation of TPR/TNR
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0,1]).ravel()

    metrics = {
        'ROC AUC' : roc_auc_score(y_true, y_proba),
        'True Negative': tn,
        'False Positive': fp,
        'False Negative': fn,
        'True Positive': tp
    }

    return metrics

def run_single_dataset(filepath, output_file, pred_col, source_col, adjustment, classifier, clf_model, n_splits, current_repeat):
    """Generate metrics for a single dataset."""
    # Load pandas dataframe
    df = pd.DataFrame()
    try:
        print_now(f"Loading combined data from {filepath}")
        df = load_dataframe(filepath, pred_col, source_col) # Give df an index - double check
    except Exception as e:
        print_now(f"Error loading or validating data: {e}")

    source = df[source_col].unique()
    if len(source) != 1:
        raise ValueError(f"Expected exactly 2 unique sources in '{source_col}', found: {source}")

    # Separate data into training and testing
    y = df[pred_col]
    cols_to_drop = [pred_col, source_col]
    X = df.drop(columns=[col for col in cols_to_drop if col in df.columns]).select_dtypes(include=[np.number])

    # Create dataframe for y predictions
    predictions = pd.DataFrame(index=df.index, columns=['y_predicted', 'y_probability'])

    # Clone the model
    model = clone(clf_model)
    cv_random_seed = 42 + current_repeat
    cv = RepeatedStratifiedKFold(n_splits=n_splits, random_state=cv_random_seed)
    splits = list(cv.split(X,y))

    # Fit the model, test, and calculate metrics
    for train_index, test_index in splits:
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        model.fit(X_train, y_train)
    
        # Predict off of x test to get y test predictions
        y_pred = model.predict(X_test)
        predictions.loc[test_index, 'y_predicted'] = y_pred

        # Take the second column of probabilistic predictions from x test
        y_proba = model.predict_proba(X_test)[:,1]
        predictions.loc[test_index, 'y_probability'] = y_proba

    # Generate the metrics from y, y predictions, and y probabilistic predictions
    metrics = calculate_metrics(y, predictions['y_predicted'], predictions['y_probability'])
    metrics_df = pd.DataFrame([metrics])

    # Add other columns
    
    metrics_df['Classifier'] = classifier
    metrics_df['Adjustment'] = adjustment
    metrics_df['Prediction'] = pred_col
    metrics_df['Train'] = source
    metrics_df['Test'] = source

    # Reorder columns and append to output file
    output_cols = ['Train', 'Test', 'ROC AUC', 'True Negative', 'False Negative', 'False Positive', 'True Positive', 
    'Classifier', 'Adjustment', 'Prediction']
    safe_write_to_csv(metrics_df[output_cols], output_file)
    print_now(f"Classification results for {filepath} saved.")

def run_combined_dataset(filepath, output_file, pred_col, source_col, adjustment, classifier, clf_model, n_splits, current_repeat):
    """Generate metrics for combined datasets
    with each training and testing combination."""
    # Load pandas dataframe
    df = pd.DataFrame()
    try:
        df = load_dataframe(filepath, pred_col, source_col)
    except Exception as e:
        print_now(f"Error loading or validating data: {e}")
        traceback.print_exc()

    # Ensure source column has two unique values
    sources = df[source_col].unique()
    if len(sources) != 2:
        raise ValueError(f"Expected exactly 2 unique sources in '{source_col}', found: {sources}")
    source1, source2 = sources

    # Define train/test combinations
    combinations = [
        ('combined', 'combined'),
        ('combined', source1),
        ('combined', source2),
        (source1, source2),
        (source2, source1)
    ]

    # Drop unwanted columns to get features
    cols_to_drop = [pred_col, source_col]
    feature_df = df.drop(columns=[col for col in cols_to_drop if col in df.columns]).select_dtypes(include=[np.number])

    for train_key, test_key in combinations:
        # Select training data
        if train_key == 'combined':
            train_df = df
        else:
            train_df = df[df[source_col] == train_key]

        # Select testing data
        if test_key == 'combined':
            test_df = df
        else:
            test_df = df[df[source_col] == test_key]

        # Get features and targets
        X_train = feature_df.loc[train_df.index]
        y_train = train_df[pred_col]

        X_test = feature_df.loc[test_df.index]
        y_test = test_df[pred_col]
    
        # Create empty dataframe for y predictions
        predictions = pd.DataFrame(index=df.index, columns=['y_predicted', 'y_probability'])

        # Clone model
        model = clone(clf_model)

        # Fit the model and create predictions
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        y_proba = model.predict_proba(X_test)[:,1]

        predictions.loc[y_test.index, 'y_predicted'] = y_pred
        predictions.loc[y_test.index, 'y_probability'] = y_proba

        y_true = y_test
        y_pred_all = y_pred
        y_proba_all = y_proba

        # Calculate metrics
        metrics = calculate_metrics(y_true, y_pred_all, y_proba_all)
        metrics_df = pd.DataFrame([metrics])
    
        # Add other columns
        metrics_df['Classifier'] = classifier
        metrics_df['Adjustment'] = adjustment
        metrics_df['Prediction'] = pred_col
        metrics_df['Train'] = train_key
        metrics_df['Test'] = test_key

        # Reorder columns and append to output file
        output_cols = ['Train', 'Test', 'ROC AUC', 'True Negative', 'False Negative', 'False Positive', 'True Positive', 
        'Classifier', 'Adjustment', 'Prediction']
        safe_write_to_csv(metrics_df[output_cols], output_file)
        print_now(f"Classification results for {filepath} saved.")

def load_dataframe(filename, pred_col, source_col):
    """Read the file into a pandas dataframe and check it has the required columns."""
    print_now(f"Loading data from {filename}")
    df = pd.read_csv(filename)

    if pred_col not in df.columns or source_col not in df.columns:
        raise ValueError(f"{pred_col} and {source_col} are required columns.")

    return df

def generate_runs(single_list, combined_list):
    """Return a list of dictionaries of run parameters for each combination of datasets."""
    
    runs = []
    for filepath in single_list:
        runs.append({
            "type": "single",
            "filename": filepath
        })

    for filepath in combined_list:
        runs.append({
            "type": "combined",
            "filename": filepath
        })

    return runs

def execute_run(args, run, model, current_repeat):
    """Perform a single run based on its type."""
    if run['type'] == 'single':
        run_single_dataset(
            filepath=run['filename'],
            output_file=args.output,
            pred_col=args.prediction_column,
            source_col=args.source_column,
            adjustment=args.adjustment,
            classifier=args.classifier,
            clf_model=model,
            n_splits=args.n_splits,
            current_repeat=current_repeat
        )

    elif run['type'] == 'combined':
        run_combined_dataset(
            filepath=run['filename'],
            output_file=args.output,
            pred_col=args.prediction_column,
            source_col=args.source_column,
            adjustment=args.adjustment,
            classifier=args.classifier,
            clf_model=model,
            n_splits=args.n_splits,
            current_repeat=current_repeat
        )

def initialize_model(classifier):
    random_seed = 42
    if classifier == 'HistGradientBoosting':
        model = HistGradientBoostingClassifier(max_iter=100, random_state=random_seed)
    elif classifier == 'RandomForest':
        model = RandomForestClassifier(n_estimators=100, n_jobs=1, random_state=random_seed)
    else:
        print_now(f"{classifier} is not currently recognized. Please try 'HistGradientBoosting' or 'RandomForest'.")

    return model

def main():
    """Parse arguments and run the classification pipeline."""
    parser = argparse.ArgumentParser(description="Run HistGradientBoosting classification on gene expression for ER status.")
    parser.add_argument('--combined-list', nargs='*', required=True, help='List of combined gene expression data filenames.')
    parser.add_argument('--single-list', nargs='*', required=True, help='List of individual gene expression data filenames.')
    parser.add_argument('--output', required=True, help='Path for the detailed output CSV file.')
    parser.add_argument('--adjustment', required=True, help='Type of adjustment done to CSV file.')
    parser.add_argument('--classifier', default= 'HistGradientBoosting', help='Name of classifier algorithm used.')
    parser.add_argument('--prediction-column', default='meta_er_status', help='Name of the prediction column, y.')
    parser.add_argument('--source-column', default='meta_source', help='Name of the source column.')
    parser.add_argument('--n', '--n-repeats', type=int, default=10, help='Number of repeats for cross-validation.')
    parser.add_argument('--n-splits', type=int, default=3, help='Number of splits for cross-validation.')
    parser.add_argument('--force-rerun', action='store_true', help='Force re-computation even if cache is valid.')
    args = parser.parse_args()

    # Create output file from output_file path and create header
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    if not os.path.exists(args.output) or os.path.getsize(args.output) == 0:
        pd.DataFrame(columns=['Train', 'Test', 'ROC AUC', 'True Negative', 'False Negative', 'False Positive', 'True Positive', 
    'Classifier', 'Adjustment', 'Prediction']).to_csv(args.output, index=False)
    
    print_now(f"=== HistGradientBoostingClassifier ER Status Classification ({args.adjustment}) ===")

    # Set up model
    model = initialize_model(args.classifier)
    
    # Generate a list of dictionaries with the parameters for each run
    runs = generate_runs(args.single_list, args.combined_list)

    # Execute runs in parallel

    n_jobs = min(len(runs), os.cpu_count() or 1)

    for current_repeat in range(args.n):
        print_now(f"=== Repeat {current_repeat + 1} of {args.n} ===")
        Parallel(n_jobs=n_jobs)(
            delayed(execute_run)(args, run, model, current_repeat) for run in runs
        )

    # for run in runs:
    #     execute_run(run)

    print_now("\nPipeline finished.")
    print_now(f"Detailed results are in: {args.output}")
    print_now("="*60)

if __name__ == "__main__":
    main()

    