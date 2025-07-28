# /scripts/metrics/reduce.py
#
# This script performs dimensionality reduction (PCA, LDA, t-SNE, or UMAP)
# on a given input data file and saves the resulting 2D coordinates and
# relevant metadata to a new CSV file. This output is intended to be used
# by a subsequent plotting script.

import pandas as pd
import numpy as np
import argparse
import os
import warnings

# Dimensionality Reduction & Scaling
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.manifold import TSNE
import umap.umap_ as umap

# Suppress warnings for a cleaner output
warnings.filterwarnings('ignore')

def run_reduction(input_file, output_file, batch_col, true_col, reduction_type="pca", debug=False):
    """
    Handles loading data, running dimensionality reduction, and saving the coordinates to a CSV.

    Args:
        input_file (str): Path to the input CSV file.
        output_file (str): Path where the output CSV will be saved.
        batch_col (str): Name of the column with batch labels.
        true_col (str): Name of the column with true biological labels.
        reduction_type (str): The reduction method ('pca', 'lda', 'tsne', 'umap').
        debug (bool): Flag to enable verbose debug messages.
    """
    if not os.path.exists(input_file):
        if debug: print(f"DEBUG: File not found, skipping: {input_file}")
        return

    if debug: print(f"DEBUG: Processing: {input_file} with {reduction_type.upper()}")

    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        if debug: print(f"DEBUG: Could not read {input_file}. Error: {e}")
        return

    # --- Data Preparation ---
    # Identify feature columns (those not starting with 'meta_')
    feature_cols = [col for col in df.columns if not col.startswith('meta_')]
    feature_data = df[feature_cols].apply(pd.to_numeric, errors='coerce')

    # Remove columns with zero or near-zero variance to avoid issues with scaling/reduction
    variances = feature_data.var()
    low_variance_cols = variances[variances < 1e-10].index
    if not low_variance_cols.empty:
        if debug: print(f"DEBUG: Removing {len(low_variance_cols)} columns with zero variance.")
        feature_data = feature_data.drop(columns=low_variance_cols)

    # Scale data before dimensionality reduction
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(feature_data)

    # --- Dimensionality Reduction ---
    dim_coords = None

    if reduction_type == 'pca':
        model = PCA(n_components=2)
        dim_coords = model.fit_transform(scaled_data)

    elif reduction_type == 'lda':
        true_labels = df[true_col]
        # LDA requires at least two classes for discrimination
        if true_labels.nunique() < 2:
            if debug: print(f"DEBUG: Skipping LDA for {input_file} - requires at least 2 distinct groups.")
            return
        model = LinearDiscriminantAnalysis(n_components=1) # n_components cannot be larger than n_classes - 1
        dim_coords = model.fit_transform(scaled_data, true_labels)

    elif reduction_type == 'tsne':
        # Use a fixed random_state for reproducibility
        model = TSNE(n_components=2, perplexity=30, random_state=42, n_iter=300)
        dim_coords = model.fit_transform(scaled_data)

    elif reduction_type == 'umap':
        # Use a fixed random_state for reproducibility
        model = umap.UMAP(n_components=2, random_state=42)
        dim_coords = model.fit_transform(scaled_data)

    if dim_coords is None:
        if debug: print(f"DEBUG: Dimensionality reduction failed for {input_file}")
        return

    # --- Create Result DataFrame and Save ---
    # The original error was a SyntaxError on the following line.
    # The fix is to provide a list of column names.
    result_df = pd.DataFrame(dim_coords)
    # Corrected column naming to handle both 1D (LDA with 2 classes) and 2D results
    result_df.columns = [f'Dim{i+1}' for i in range(result_df.shape[1])] # Changed column assignment. Might fix column mismatch error.

    # If LDA produced only one dimension, duplicate it for plotting purposes
    if result_df.shape[1] == 1:
        result_df['Dim2'] = result_df['Dim1']


    # Add the metadata columns needed for plotting later.
    result_df[batch_col] = df[batch_col].values
    result_df[true_col] = df[true_col].values

    # Save to output file
    try:
        # Create the directory for the output file if it doesn't exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        result_df.to_csv(output_file, index=False)
        if debug: print(f"DEBUG: Successfully saved reduced data to {output_file}")
    except Exception as e:
        if debug: print(f"DEBUG: Could not save to {output_file}. Error: {e}")


def main():
    """ Main execution block for command-line use """
    parser = argparse.ArgumentParser(description="Run dimensionality reduction and save coordinates to a CSV file.")
    parser.add_argument("--input-file", type=str, required=True, help="Input CSV file with expression data.")
    parser.add_argument("--output-file", type=str, required=True, help="Path to save the output CSV with reduced coordinates.")
    parser.add_argument("--batch-col", type=str, required=True, help="Name of the column with batch information.")
    parser.add_argument("--true-col", type=str, required=True, help="Name of the column with the true biological signal.")
    parser.add_argument("--reduction-type", type=str, default="pca", choices=["pca", "lda", "tsne", "umap"], help="Type of reduction to perform.")
    parser.add_argument("--debug", action="store_true", help="Enable debug messages.")
    args = parser.parse_args()

    run_reduction(
        input_file=args.input_file,
        output_file=args.output_file,
        batch_col=args.batch_col,
        true_col=args.true_col,
        reduction_type=args.reduction_type,
        debug=args.debug
    )

if __name__ == "__main__":
    main()
