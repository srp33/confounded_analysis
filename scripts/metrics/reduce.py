# /scripts/metrics/reduce.py
#
# This script performs all dimensionality reductions (PCA, LDA, t-SNE, UMAP)
# on a given input data file. It uses a granular, hash-based caching system
# to avoid re-processing individual reduction steps that have not changed.

import pandas as pd
import numpy as np
import argparse
import os
import warnings
from pathlib import Path

# --- Local Imports ---
import sys
import os
# Add the parent directory (scripts) to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from scripts.utils import HashCache
except ImportError:
    from utils import HashCache

# Suppress TensorFlow informational messages
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
warnings.filterwarnings('ignore')

# Dimensionality Reduction & Scaling
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.manifold import TSNE
import umap.umap_ as umap
from sklearn.impute import SimpleImputer


def print_now(*args, **kwargs):
    print(*args, **kwargs, flush=True)


def run_all_reductions(hash_cache, input_file, output_dir, batch_col, true_col, additional_meta_cols=None, debug=False):
    """
    Handles loading data once, then running and caching each dimensionality
    reduction method individually.
    """
    if not os.path.exists(input_file):
        if debug: print_now(f"DEBUG: File not found, skipping: {input_file}")
        return

    method_name = os.path.splitext(os.path.basename(input_file))[0]
    print_now(f"  - Processing file: {input_file} ")

    # --- Data Preparation (Done Once) ---
    # This block is executed only once per input file, before the loop.
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        if debug: print_now(f"DEBUG: Could not read {input_file}. Error: {e}")
        return

    required_cols = [batch_col, true_col]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        if debug: print_now(f"DEBUG: ERROR - Required metadata columns not found: {missing_cols} in {input_file}")
        return
    
    # Prepare list of all metadata columns to include
    meta_cols_to_include = [batch_col, true_col]
    if additional_meta_cols:
        for col in additional_meta_cols:
            if col in df.columns and col not in meta_cols_to_include:
                meta_cols_to_include.append(col)

    feature_cols = [col for col in df.columns if not col.startswith('meta_')]
    feature_data = df[feature_cols].apply(pd.to_numeric, errors='coerce')

    if feature_data.isnull().values.any():
        if debug: print_now(f"DEBUG: Input contains NaN values. Imputing with column means. File: {input_file}")
        imputer = SimpleImputer(strategy='mean', keep_empty_features=True)
        imputed_data = imputer.fit_transform(feature_data)
        feature_data = pd.DataFrame(imputed_data, index=feature_data.index, columns=feature_data.columns)

    variances = feature_data.var()
    low_variance_cols = variances[variances < 1e-10].index
    if not low_variance_cols.empty:
        if debug: print_now(f"DEBUG: Removing {len(low_variance_cols)} columns with zero/low variance.")
        feature_data = feature_data.drop(columns=low_variance_cols)
    
    if feature_data.shape[1] == 0:
        if debug: print_now(f"DEBUG: No feature columns with variance remain. Skipping {input_file}.")
        return

    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(feature_data)

    # --- Loop Through and Cache Each Reduction Type Individually ---
    # The 'lda' reduction now uses the combined batch and true class labels.
    reduction_types = ["pca", "lda", "tsne", "umap"]
    for reduction_type in reduction_types:
        cache_key = f"reduction|{method_name}|{reduction_type}"
        output_file = os.path.join(output_dir, f"{method_name}-{reduction_type}.csv")

        # The context manager handles hash checking and updating.
        # The block is only entered if the file needs processing.
        with hash_cache.check(cache_key, input_file) as should_skip:
            # Also check if the output file exists, in case it was manually deleted.
            if should_skip and os.path.exists(output_file):
                if debug: print_now(f"    Skipping {reduction_type.upper()} for {input_file} (no changes)")
                continue

            # This block runs only if the file is new, has changed, or the output is missing.
            if debug: print_now(f"    -- Running {reduction_type.upper()} for {input_file}")
            dim_coords = None
            try:
                if reduction_type == 'pca':
                    model = PCA(n_components=2)
                    dim_coords = model.fit_transform(scaled_data)
                elif reduction_type == 'lda':
                    # --- LDA Grouping Logic ---
                    # Combine batch and true labels to form distinct groups for LDA.
                    # e.g., ('Batch1', 'True') -> 'Batch1_True'
                    combined_labels = df[batch_col].astype(str) + '_' + df[true_col].astype(str)
                    
                    # Ensure there are at least two unique groups to compare.
                    if combined_labels.nunique() < 2:
                        if debug: print_now(f"DEBUG: Skipping LDA - requires at least 2 distinct combined groups.")
                        continue
                    
                    # LDA components are limited by min(n_features, n_classes - 1).
                    # We want 2 components, capped by the number of available groups.
                    n_components = min(2, combined_labels.nunique() - 1)

                    model = LinearDiscriminantAnalysis(n_components=n_components)
                    dim_coords = model.fit_transform(scaled_data, combined_labels)
                elif reduction_type == 'tsne':
                    model = TSNE(n_components=2, perplexity=min(30, len(df)-1), random_state=42, max_iter=300)
                    dim_coords = model.fit_transform(scaled_data)
                elif reduction_type == 'umap':
                    model = umap.UMAP(n_components=2, random_state=42)
                    dim_coords = model.fit_transform(scaled_data)
            except Exception as e:
                if debug: print_now(f"DEBUG: {reduction_type.upper()} failed for {input_file}. Error: {e}")
                continue # Skip to the next reduction type

            # --- Create Result DataFrame and Save ---
            if dim_coords is None:
                continue # Skip if reduction failed or was skipped

            result_df = pd.DataFrame(dim_coords)
            result_df.columns = [f'Dim{i+1}' for i in range(result_df.shape[1])]
            
            # If LDA produces only one dimension, create a second dimension of zeros
            # to maintain a consistent 2D output format.
            if result_df.shape[1] == 1:
                result_df['Dim2'] = 0
            
            # Add all metadata columns
            for col in meta_cols_to_include:
                result_df[col] = df[col].values
            
            try:
                os.makedirs(output_dir, exist_ok=True)
                result_df.to_csv(output_file, index=False)
            except Exception as e:
                if debug: print_now(f"DEBUG: Could not save to {output_file}. Error: {e}")


def main():
    """ Main execution block for command-line use """
    parser = argparse.ArgumentParser(description="Run all dimensionality reductions on an input file and save coordinate CSVs.")
    parser.add_argument("--input-file", type=str, required=True, help="Input CSV file with expression data.")
    parser.add_argument("--output-dir", type=str, required=True, help="Path to save the output CSVs with reduced coordinates.")
    parser.add_argument("--batch-col", type=str, required=True, help="Name of the column with batch information.")
    parser.add_argument("--true-col", type=str, required=True, help="Name of the column with the true biological signal.")
    parser.add_argument("--additional-meta-cols", type=str, nargs='*', help="Additional metadata columns to include in output.")
    parser.add_argument("--hash-dir", type=Path, required=True, help="Directory to store cache files.")
    parser.add_argument("--write-over", action="store_true", help="Force re-computation and overwrite cache.")
    parser.add_argument("--debug", action="store_true", help="Enable detailed debug prints.")
    args = parser.parse_args()

    # Determine the dataset name from the input file path
    # e.g., /data/gold/gse20194/combat.csv -> gse20194
    try:
        dataset_name = Path(args.input_file).parts[-2]
    except IndexError:
        dataset_name = "unknown_dataset"

    hash_cache = HashCache(
        hash_dir=args.hash_dir,
        cache_filename=f"{dataset_name}_reduce.hashes.json",
        write_over=args.write_over,
        debug=args.debug
    )

    run_all_reductions(
        hash_cache=hash_cache,
        input_file=args.input_file,
        output_dir=args.output_dir,
        batch_col=args.batch_col,
        true_col=args.true_col,
        additional_meta_cols=args.additional_meta_cols,
        debug=args.debug
    )

if __name__ == "__main__":
    main()