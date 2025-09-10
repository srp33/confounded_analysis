import argparse
import os
import sys
import pandas as pd
from pathlib import Path
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, confusion_matrix
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.base import clone
from joblib import Parallel, delayed

# Add the parent directory (scripts) to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def print_now(*args, **kwargs):
    """Print a message to the console with flushing."""
    print(*args, flush=True, **kwargs)

def run_single_dataset(filename, pred_col, source_col, model, n_splits, n_repeats, current_repeat):
    """Generate metrics for a single dataset."""
    # Load pandas dataframe
    try:
        df = load_dataframe(filename, pred_col, source_col) # Give df an index - double check
    except Exception as e:
        print_now(f"Error loading or validating data: {e}")

    if n_repeats <= 0:
        return
    print_now(f"\nRunning classification for {filename} with {n_repeats} repeats (offset: {current_repeat})...")

    # Separate data into training and testing
    y = df[pred_col]
    cols_to_drop = [pred_col, source_col]
    X = df.drop(columns=[col for col in cols_to+drop if col in df.columns]).select_dtypes(include=[np.number])

    # Clone the model and create random seed for the cross-validation with repeat offset
    random_seed = 42 + repeat_offset
    model = HistGradientBoostingClassifier(max_iter=100, random_state=random_seed)
    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=random_seed)
    splits = list(cv.split(X,y))

    # Create empty dataframe - match the index of y (before train or test) or create another column
    new_df = pandas.DataFrame(index=df.index)

    for train_index, test_index in splits:
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        model.fit(X_train, y_train)
    
        # Predict off of x test to get y test predictions
        # Probabilistic predictions from x test
        new_d.ilocf[test_index] = probabilistic_predictions

    # Generate the metrics from y, y predictions, and y probabilistic predictions

    # Create a dataframe of the metrics

def run_combined_dataset(filename, pred_col, source_col, model, n_splits, n_repeats, current_repeat):
    """Generate metrics for combined datasets
    with each training and testing combination."""
    # Load pandas dataframe
    try:
        df = load_dataframe(filename, pred_col, source_col)
    except Exception as e:
        print_now(f"Error loading or validating data: {e}")

    if n_repeats <= 0:
        return
    print_now(f"\nRunning classification for {filename} with {n_repeats} repeats (offset: {current_repeat})...")


    # Separate data into training and testing

    # Clone the model and fit it to the data

    # Create a dataframe of the metrics

def load_dataframe(filename, pred_col, source_col):
    """Read the file into a pandas dataframe and check it has the required columns."""
    print_now(f"Loading data from {filename}")
    df = pd.read_csv(filename)

    if pred_col not in df.columns or source_col not in df.columns:
        raise ValueError(f"{pred_col} and {source_col} are required columns.")

    return df

def main():
    """Parse arguments and run the classification pipeline."""
    parser = argparse.ArgumentParser(description="Run HistGradientBoosting classification on gene expression for ER status.")
    parser.add_argument('--combined-list', nargs='*', required=True, help='List of combined gene expression data filenames.')
    parser.add_argument('--single-list', nargs='*', required=True, help='List of individual gene expression data filenames.')
    parser.add_argument('--output', required=True, help='Path for the detailed output CSV file.')
    parser.add_argument('--confusion-matrix', required=True, help='Path for the confusion matrix .txt file.')
    parser.add_argument('--summary', required=True, help='Path for the summary metrics CSV file.')
    parser.add_argument('--adjustment', default='Unadjusted', help='Name of the adjustment method used.')
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
        pd.DataFrame(columns=['train', 'test', 'ROC_AUC', 'true_negative', 'false_positive', 'false_negative', 'true positive',
        'classifier', 'adjuster', 'dataset', 'pred_column'])
    
    print_now(f"=== HistGradientBoostingClassifier ER Status Classification ({args.adjustment}) ===")

    # Set up model and cross-validation
        # Set a random seed
        # Define the model (can be variable)

    # Loop to evaluate each dataset combination, appending its returned metrics dataframe

    # Aggregate results and add to .csv file