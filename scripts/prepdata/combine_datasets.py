# combine_datasets.py - Simplified version with gene ID conversion moved upstream

import argparse
import os
import pandas as pd
from sklearn.preprocessing import LabelEncoder

def print_now(*args, **kwargs):
    """Prints a message to the console with flushing to ensure immediate output."""
    print(*args, flush=True, **kwargs)

def load_and_prepare_data(file_path, debug=False):
    """
    Loads a single dataset and prepares it for combination.
    
    Note: Gene ID conversion is now handled upstream in convert_raw_files.py,
    so datasets should already have gene symbols as column names.

    Args:
        file_path (str): Path to the input CSV file.
        debug (bool): If True, enables detailed print statements for debugging.

    Returns:
        pd.DataFrame: The processed dataframe.
    """
    print_now(f"Loading data from {file_path}")
    df = pd.read_csv(file_path, low_memory=False)

    # Ensure the target column 'meta_er_status' exists.
    if 'meta_er_status' not in df.columns:
        raise ValueError(f"'meta_er_status' column not found in {file_path}")

    if debug:
        print_now(f"DEBUG: Loaded dataset with shape: {df.shape}")
        gene_cols = [col for col in df.columns if not col.startswith('meta_')]
        print_now(f"DEBUG: Found {len(gene_cols)} gene columns")
        print_now(f"DEBUG: Example gene columns: {gene_cols[:5]}")

    return df

def get_common_genes(df1, df2, debug=False):
    """
    Finds the intersection of gene columns between two dataframes.

    Args:
        df1 (pd.DataFrame): The first dataframe.
        df2 (pd.DataFrame): The second dataframe.
        debug (bool): If True, enables detailed print statements for debugging.

    Returns:
        list: A list of common gene symbols.
    """
    # Extract gene columns, excluding metadata and sample identifiers.
    genes1 = [col for col in df1.columns if not col.startswith('meta_') and col != "Sample_ID"]
    genes2 = [col for col in df2.columns if not col.startswith('meta_') and col != "Sample_ID"]
    
    # Find the intersection of the two gene lists.
    common_genes = list(set(genes1) & set(genes2))
    print_now(f"Found {len(common_genes)} common genes between datasets.")
    if debug and common_genes:
        print_now(f"DEBUG: Example common genes: {common_genes[:5]}")
    
    if not common_genes:
        print_now("ERROR: No common genes found between datasets. Please check your input files.")
        print_now(f"First 10 genes of dataset 1: {genes1[:10]}")
        print_now(f"First 10 genes of dataset 2: {genes2[:10]}")
    elif len(common_genes) < 10:
        print_now(f"WARNING: Very few common genes found. Common genes: {common_genes}")
    return common_genes

def main():
    """Main function to parse arguments and orchestrate the dataset combination."""
    parser = argparse.ArgumentParser(
        description="Combine two gene expression datasets based on common genes.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    # --- Input/Output Arguments ---
    parser.add_argument('--input1', required=True, help='Path to the first data file.')
    parser.add_argument('--input2', required=True, help='Path to the second data file.')
    parser.add_argument('--output', required=True, help='Path for the combined output CSV file.')

    # --- General Arguments ---
    parser.add_argument('--debug', action='store_true', help='Enable detailed debug output.')

    args = parser.parse_args()

    try:
        # --- Load and Prepare Datasets ---
        print_now("--- Processing Dataset 1 ---")
        df1 = load_and_prepare_data(file_path=args.input1, debug=args.debug)
        
        print_now("\n--- Processing Dataset 2 ---")
        df2 = load_and_prepare_data(file_path=args.input2, debug=args.debug)
    except Exception as e:
        print_now(f"Error during data loading or preparation: {e}")
        return

    # --- Combine Datasets ---
    print_now("\n--- Combining Datasets ---")
    common_genes = get_common_genes(df1, df2, debug=args.debug)
    if not common_genes:
        print_now("No common genes found. Cannot create a combined dataset.")
        return

    # Add a source identifier to each dataset using the folder name
    source1_name = os.path.basename(os.path.dirname(args.input1))
    source2_name = os.path.basename(os.path.dirname(args.input2))
    print_now(f"Source names: {source1_name}, {source2_name}")

    df1['meta_source'] = source1_name
    df2['meta_source'] = source2_name

    # Define columns to keep for the final combined dataframe.
    keep_cols = common_genes + ['meta_er_status', 'meta_source']
    
    # Filter each dataframe to common columns and concatenate them.
    df_combined = pd.concat(
        [df1[keep_cols], df2[keep_cols]], 
        ignore_index=True
    )
    
    print_now(f"Combined dataset shape: {df_combined.shape}")

    # --- Save Combined Dataset ---
    print_now(f"\n--- Saving Combined Dataset ---")
    
    # Ensure output directory exists
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    df_combined.to_csv(args.output, index=False)
    print_now(f"Combined dataset saved to: {args.output}")
    
    # Final summary
    print_now(f"\n--- Summary ---")
    print_now(f"Dataset 1: {df1.shape[0]} samples")
    print_now(f"Dataset 2: {df2.shape[0]} samples")
    print_now(f"Combined: {df_combined.shape[0]} samples, {len(common_genes)} genes")
    print_now(f"Output: {args.output}")

if __name__ == "__main__":
    main()