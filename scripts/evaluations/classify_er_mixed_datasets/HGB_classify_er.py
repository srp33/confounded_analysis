import argparse
import os
import sys
import pandas as pd

def print_now(*args, **kwargs):
    """Print a message to the console with flushing."""
    print(*args, flush=True, **kwargs)

def run_single_dataset(filename, pred_col, source_col, model, current_repeat_num, num_splits):
    """Generate metrics for a single dataset."""
    # Create pandas dataframe

    # Separate data into training and testing

    # Clone the model and fit it to the data

    # Create a dataframe of the metrics

def run_combined_dataset(filename, pred_col, source_col, model, current_repeat_num, num_splits):
    """Generate metrics for combined datasets
    with each training and testing combination."""
    # Create pandas dataframe

    # Separate data into training and testing

    # Clone the model and fit it to the data

    # Create a dataframe of the metrics

def generate_summary(output_file, summary_file):
    """Read detailed results and generate a summary with mean, std, and SEM."""
    print_now("\n" + "="*60 + "\nGENERATING SUMMARY METRICS")
    try:
        df = pd.read_csv(output_file)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        print_now(f"Warning: Detailed results file '{output_file}' not found or empty. Cannot generate summary.")
        return

    # Group by experiment identifiers and calculate aggregate stats.
    summary = df.groupby(['adjuster', 'dataset', 'metric'])['value'].agg(['mean', 'std']).reset_index()
    summary['sem'] = df.groupby(['adjuster', 'dataset', 'metric'])['value'].apply(lambda x: x.std() / np.sqrt(x.count()) if x.count() > 1 else 0.0).values
    
    summary.rename(columns={'dataset': 'evaluation'}, inplace=True)
    
    os.makedirs(os.path.dirname(summary_file), exist_ok=True)
    summary.to_csv(summary_file, index=False, float_format='%.4f')
    print_now(f"Summary metrics saved to: {summary_file}")
    print_now("="*60)

def print_confusion_matrix(output_file, matrix_file):
    print_now("Printing confusion matrix")
    # Filter confusion matrix values from metrics CSV file
    metrics_full = pd.read_csv(output_file)
    matrix_values = ['True Negative', 'False Positive', 'False Negative', 'True Positive']
    metrics_filtered = metrics_full[metrics_full['combination'].isin(matrix_values)]
    
    # Create the .txt file
    os.makedirs(os.path.dirname(matrix_file), exist_ok=True)
    metrics_filtered.to_csv(matrix_file, sep='\t', index=False)
    print_now(f"Confusion matrix values safed to: {matrix_file}")

def main(dataset_list, output_file, classifier, pred_col, source_col,  num_splits, num_folds):
    """Parse arguments and run the classification pipeline."""
    # Define arguments

    # Create output file from output_file path and create header

    # Create a list of single and a list of combined dataset combinations from the provided dataset_list

    # Set up model and cross-validation

    # Loop to evaluate each dataset combination, appending its returned metrics dataframe

    # Aggregate results and add to .csv file