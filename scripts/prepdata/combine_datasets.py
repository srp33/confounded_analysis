# combine_datasets.py

import argparse
import os
import pandas as pd
from sklearn.preprocessing import LabelEncoder

def print_now(*args, **kwargs):
    """Prints a message to the console with flushing to ensure immediate output."""
    print(*args, flush=True, **kwargs)

def map_probes_to_genes(df, annotation_file, debug=False):
    """
    Maps Affymetrix probe IDs to gene symbols and aggregates expression by gene.

    Args:
        df (pd.DataFrame): The input dataframe with probe IDs as columns.
        annotation_file (str): Path to the CSV annotation file mapping ProbeIDs to GeneSymbols.
        debug (bool): If True, enables detailed print statements for debugging.

    Returns:
        pd.DataFrame: A new dataframe with gene symbols as columns.
    """
    print_now(f"Mapping probe IDs to gene symbols using {annotation_file}...")

    # Load annotation, standardizing probe IDs by removing common suffixes.
    annot = pd.read_csv(annotation_file)
    annot['ProbeID'] = annot['ProbeID'].str.replace(r'_[a-z]_at$', '_at', regex=True)

    if debug:
        print_now(f"DEBUG: Annotation file loaded with shape: {annot.shape}")

    # Clean annotation data and create a mapping dictionary.
    annot = annot.dropna(subset=['ProbeID', 'GeneSymbol'])
    annot = annot.drop_duplicates(subset=['ProbeID'], keep='first')
    probe_to_gene_map = annot.set_index('ProbeID')['GeneSymbol'].to_dict()

    # Separate metadata (columns starting with 'meta_') from expression data.
    meta_cols = [col for col in df.columns if col.startswith('meta_')]
    meta_df = df[meta_cols]
    probe_cols = [col for col in df.columns if not col.startswith('meta_')]
    probe_df = df[probe_cols]

    # Transpose, map probe IDs to gene symbols, and aggregate by taking the mean.
    probe_df_T = probe_df.T
    probe_df_T['GeneSymbol'] = probe_df_T.index.map(probe_to_gene_map)
    probe_df_T = probe_df_T.dropna(subset=['GeneSymbol'])
    gene_df_T = probe_df_T.groupby('GeneSymbol').mean(numeric_only=True)

    # Transpose back and merge with original metadata.
    gene_df = gene_df_T.T
    final_df = pd.concat([meta_df.reset_index(drop=True), gene_df.reset_index(drop=True)], axis=1)

    print_now(f"Probe mapping complete. New data shape: {final_df.shape}")
    return final_df

def map_entrez_to_genesymbols(df, annotation_file, debug=False):
    """
    Maps Entrez Gene IDs to gene symbols and aggregates expression by gene.

    Args:
        df (pd.DataFrame): The input dataframe with Entrez IDs as columns.
        annotation_file (str): Path to the CSV annotation file mapping EntrezID to GeneSymbol.
        debug (bool): If True, enables detailed print statements for debugging.

    Returns:
        pd.DataFrame: A new dataframe with gene symbols as columns.
    """
    print_now(f"Mapping Entrez IDs to gene symbols using {annotation_file}...")
    annot = pd.read_csv(annotation_file)

    # Clean annotation data and create a mapping dictionary.
    annot = annot.dropna(subset=['EntrezID', 'GeneSymbol'])
    entrez_to_gene_map = annot.set_index('EntrezID')['GeneSymbol'].to_dict()

    # Separate metadata from expression data.
    meta_cols = [col for col in df.columns if col.startswith('meta_')]
    meta_df = df[meta_cols]
    entrez_cols = [col for col in df.columns if not col.startswith('meta_')]
    entrez_df = df[entrez_cols]
    
    # Standardize column names by removing suffixes.
    entrez_df.columns = entrez_df.columns.str.replace('_at$', '', regex=True)

    # Transpose, map Entrez IDs to gene symbols, and aggregate by taking the mean.
    entrez_df_T = entrez_df.T
    entrez_df_T.index = entrez_df_T.index.astype(int)
    entrez_df_T['GeneSymbol'] = entrez_df_T.index.map(entrez_to_gene_map)
    entrez_df_T = entrez_df_T.dropna(subset=['GeneSymbol'])
    gene_df_T = entrez_df_T.groupby('GeneSymbol').mean(numeric_only=True)
    gene_df = gene_df_T.T
    
    # Merge with original metadata.
    final_df = pd.concat([meta_df.reset_index(drop=True), gene_df.reset_index(drop=True)], axis=1)

    print_now(f"Entrez mapping complete. New data shape: {final_df.shape}")
    return final_df

def load_and_prepare_data(file_path, map_from='none', annotation_file=None, debug=False):
    """
    Loads a single dataset, performs optional feature mapping, and prepares the target variable.

    Args:
        file_path (str): Path to the input CSV file.
        map_from (str): The type of feature ID to map from ('probe', 'entrez', or 'none').
        annotation_file (str, optional): Path to the annotation file. Required if mapping.
        debug (bool): If True, enables detailed print statements for debugging.

    Returns:
        pd.DataFrame: The processed dataframe.
    """
    print_now(f"Loading data from {file_path}")
    df = pd.read_csv(file_path)

    # Map features to gene symbols if requested.
    if map_from == 'probe' and annotation_file:
        df = map_probes_to_genes(df, annotation_file, debug=debug)
    elif map_from == 'entrez' and annotation_file:
        df = map_entrez_to_genesymbols(df, annotation_file, debug=debug)

    # Ensure the target column 'meta_er_status' exists.
    if 'meta_er_status' not in df.columns:
        raise ValueError(f"'meta_er_status' column not found in {file_path}")

    # Clean and encode the target variable.
    df = df.dropna(subset=['meta_er_status'])
    df['meta_er_status'] = LabelEncoder().fit_transform(df['meta_er_status'])
    df = df.dropna(subset=['meta_er_status'])

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

    # --- Mapping Arguments for Dataset 1 ---
    parser.add_argument('--annot1', help='Optional: Path to the annotation file for the first dataset.')
    parser.add_argument('--map_type1', choices=['probe', 'entrez'], help="Mapping type for the first dataset. Required if --annot1 is provided.")

    # --- Mapping Arguments for Dataset 2 ---
    parser.add_argument('--annot2', help='Optional: Path to the annotation file for the second dataset.')
    parser.add_argument('--map_type2', choices=['probe', 'entrez'], help="Mapping type for the second dataset. Required if --annot2 is provided.")

    # --- General Arguments ---
    parser.add_argument('--debug', action='store_true', help='Enable detailed debug printing.')
    
    args = parser.parse_args()

    # --- Argument Validation ---
    if args.annot1 and not args.map_type1:
        parser.error("--map_type1 is required when --annot1 is provided.")
    if args.annot2 and not args.map_type2:
        parser.error("--map_type2 is required when --annot2 is provided.")

    try:
        # --- Load and Prepare Datasets ---
        print_now("--- Processing Dataset 1 ---")
        df1 = load_and_prepare_data(
            file_path=args.input1, 
            map_from=args.map_type1 if args.annot1 else 'none', 
            annotation_file=args.annot1, 
            debug=args.debug
        )
        
        print_now("\n--- Processing Dataset 2 ---")
        df2 = load_and_prepare_data(
            file_path=args.input2, 
            map_from=args.map_type2 if args.annot2 else 'none', 
            annotation_file=args.annot2, 
            debug=args.debug
        )
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
    source1_name =  os.path.basename(os.path.dirname(args.input1))
    source2_name =  os.path.basename(os.path.dirname(args.input2))
    print_now(f"Source names:  {source1_name}, {source2_name}")

    df1['meta_source'] = source1_name
    df2['meta_source'] = source2_name

    # Define columns to keep for the final combined dataframe.
    keep_cols = common_genes + ['meta_er_status', 'meta_source']
    
    # Filter each dataframe to common columns and concatenate them.
    df_combined = pd.concat(
        [df1[keep_cols], df2[keep_cols]], 
        ignore_index=True
    )
    
    # --- Save Output ---
    output_dir = os.path.dirname(args.output)
    if output_dir: # Create directory only if a path is specified.
        os.makedirs(output_dir, exist_ok=True)
        
    df_combined.to_csv(args.output, index=False)
    print_now(f"\nSuccessfully created combined dataset with {len(common_genes)} genes.")
    print_now(f"Final shape: {df_combined.shape}")
    print_now(f"Saved combined data to: {args.output}")

if __name__ == "__main__":
    main()
