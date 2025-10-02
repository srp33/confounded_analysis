import pandas as pd
from pathlib import Path
import argparse
from functools import reduce
from sklearn.preprocessing import LabelEncoder

# Store a single combined dataset in data/all_combined_data/

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
    df = pd.read_csv(file_path, sep='\t')

    if map_from == 'probe' and annotation_file:
        df = map_probes_to_genes(df, annotation_file, debug)
    elif map_from == 'entrez' and annotation_file:
        df = map_entrez_to_genesymbols(df, annotation_file, debug)

    if 'meta_er_status' not in df.columns:
        raise ValueError(f"'meta_er_status' not found in {file_path.name}")
    
    df = df.dropna(subset=['meta_er_status'])
    df['meta_er_status'] = LabelEncoder().fit_transform(df['meta_er_status'])

    return df


def get_common_genes_multi(dfs):
    gene_sets = [
        set(c for c in df.columns if not c.startswith('meta_') and c!= "Sample_ID")
        for df in dfs
    ]
    return list(reduce(set.intersection, gene_sets))

def main():
    parser = argparse.ArgumentParser(description="Combine GEO expression datasets into one matrix.")
    parser.add_argument('--input-dir', type=Path, required=True, help='Directory containing expression_*.tsv files')
    parser.add_argument('--output-file', type=Path, required=True, help='Path to save combined expression file')
    parser.add_argument('--map-types', nargs='*', choices=['probe', 'entrez', 'none'], help='Mapping type per file')
    parser.add_argument('--annot-dir', type=Path, help='Directory containing annotation files')
    parser.add_argument('--debug', action='store_true')

    args = parser.parse_args()

    expression_files = sorted(args.input_dir.rglob("expression_*.tsv"))
    if not expression_files:
        raise FileNotFoundError("No expression files found.")
    
    dfs = []
    for idx, file in enumerate(expression_files):
        map_type = args.map_types[idx] if idx < len(args.map_types) else 'none'
        annot_file = None
        if map_type != 'none' and args.annot_dir:
            annot_file = args.annot_dir / f"{file.stem}_annotation.csv"

        print(f"Processing {file.name} (map: {map_type})")
        df = load_and_prepare_data(file, map_from=map_type, annotation_file=annot_file, debug=args.debug)

        df['meta_source'] = file.stem
        dfs.append(df)

    common_genes = get_common_genes_multi(dfs)
    print(f"Found {len(common_genes)} common genes across all datasets.")

    # Keep only common genes and key metadata
    keep_cols = common_genes + ['meta_er_status', 'meta_source']
    combined_df = pd.concat([df[keep_cols] for df in dfs], ignore_index=True)

    # Make sure the output folder exists
    output_dir = args.output_file.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Saving combined dataset with shape {combined_df.shape}")
    combined_df.to_csv(args.output_file, index=False, sep='\t')

if __name__ == "__main__":
    main()