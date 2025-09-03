import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from functools import reduce

# --- Global Configuration ---

# Configure pandas display options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

sns.set_theme(style="whitegrid", palette="Set2")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans', 'Bitstream Vera Sans', 'sans-serif']


# --- Function Definitions ---

def load_datasets(data_dir):
    """Load all datasets from .csv files in a specified directory."""
    datasets = {}
    print(f"Loading datasets from: {data_dir}")
    try:
        for filename in os.listdir(data_dir):
            if filename.endswith('.csv'):
                if filename == 'unadjusted_t.csv':
                    continue
                # Generate a clean name from the filename
                dataset_name = os.path.splitext(filename)[0]
                filepath = os.path.join(data_dir, filename)
                datasets[dataset_name] = pd.read_csv(filepath)
                print(f"  - Loaded '{filename}' as '{dataset_name}'")
        if not datasets:
            print("Warning: No CSV files found in the directory.")
            return None
        print("All datasets loaded successfully.")
        return datasets
    except FileNotFoundError as e:
        print(f"Error: Directory not found at {data_dir}. {e}")
        return None
    except Exception as e:
        print(f"An error occurred while loading data: {e}")
        return None

def analyze_dataframe(df, df_name):
    """Analyze a single dataframe and print detailed statistics."""
    target_genes = ['ESR1', 'GP2', 'PTK6', 'STX17', 'PPP1R3C', 'SASH1']
    all_gene_columns = [col for col in df.columns if not col.startswith('meta_')]
    gene_columns = [gene for gene in target_genes if gene in all_gene_columns]
    
    # Handle cases where target genes are not found
    missing_genes = set(target_genes) - set(all_gene_columns)
    if missing_genes:
        print(f"NOTE: The following target genes are not in {df_name}: {list(missing_genes)}")
    if not gene_columns:
        print(f"WARNING: No target genes found in {df_name}, using first 6 genes instead.")
        gene_columns = all_gene_columns[:6]
    
    print(f"\n{'='*100}")
    print(f"STATISTICAL ANALYSIS FOR: {df_name.upper()}")
    print(f"{'='*100}")
    print(f"Analyzing {len(gene_columns)} gene columns: {gene_columns}")
    print(f"Dataset shape: {df.shape}")
    print(f"Meta_source groups: {df['meta_source'].unique()}")
    
    # Grouped statistics
    print("\n" + "="*80)
    print("STATISTICS GROUPED BY META_SOURCE")
    print("="*80)
    grouped = df.groupby('meta_source')[gene_columns]
    print("\nMinimum by group:\n", grouped.min().to_string())
    print("\nMean by group:\n", grouped.mean().to_string())
    
    return gene_columns

def prepare_data_for_plotting(df, df_name, debug=False):
    """Prepare a single dataframe for plotting by adding dataset info and mapping ER status."""
    df_copy = df.copy()
    df_copy['dataset'] = df_name
    if 'meta_er_status' in df_copy.columns:
        if debug:
            print(f"DEBUG: In prepare_data_for_plotting for {df_name}:")
            print(f"  Unique 'meta_er_status' values BEFORE mapping: {df_copy['meta_er_status'].unique()}")
        
        # Map numeric ER status to string representation
        df_copy['er_status_str'] = df_copy['meta_er_status'].map({1: 'Positive', 0: 'Negative'})
        df_copy.dropna(subset=['er_status_str'], inplace=True)
        
        if debug:
            print(f"  Value counts for 'er_status_str' AFTER mapping:\n{df_copy['er_status_str'].value_counts(dropna=False)}")

    elif debug:
        print(f"DEBUG: 'meta_er_status' column NOT FOUND in {df_name}")

    return df_copy

def generate_histograms(datasets, output_dir="/outputs/figures/histograms/", debug=False):
    """Generate log-transformed KDE plots to compare gene distributions."""
    print(f"\n{'='*100}")
    print("GENERATING DISTRIBUTION PLOTS (KDE)")
    print(f"{'='*100}")
    os.makedirs(output_dir, exist_ok=True)

    target_genes = ['ESR1', 'GP2', 'PTK6', 'STX17', 'PPP1R3C', 'SASH1']
    
    # Find common genes across all datasets
    common_genes = list(reduce(set.intersection, [set(df.columns) for df in datasets.values()]))
    genes_to_plot = [gene for gene in target_genes if gene in common_genes]
    
    print(f"Generating plots for {len(genes_to_plot)} common target genes: {genes_to_plot}")

    # --- FIX STARTS HERE ---
    # # Pending, might fix blank plot error (The original preparation function removed rows needed by histograms).
    # # The old method filtered data based on 'er_status_str', which is not used here and could remove all data for a facet.
    # # This new approach performs a simpler preparation, only adding the 'dataset' column.
    prepared_dfs = []
    for name, df in datasets.items():
        df_copy = df.copy()
        df_copy['dataset'] = name
        prepared_dfs.append(df_copy)
    combined_df = pd.concat(prepared_dfs, ignore_index=True)
    
    if debug:
        print(f"DEBUG: Shape of combined_df for histograms: {combined_df.shape}")
        print(f"DEBUG: Value counts for 'dataset' column:\n{combined_df['dataset'].value_counts(dropna=False)}")
    # --- FIX ENDS HERE ---
    
    for gene in genes_to_plot:
        print(f"  - Processing {gene}...")
        
        # Apply log transform safely to a temporary column
        log_gene_col = f"log_{gene}"
        # Ensure the column has no NaNs before transformation to avoid warnings/errors
        valid_data = combined_df[[gene]].dropna()
        combined_df[log_gene_col] = np.log1p(combined_df[gene] - valid_data[gene].min())
        
        # Create a faceted KDE plot for better comparison
        g = sns.displot(
            data=combined_df, 
            x=log_gene_col, 
            hue='meta_source', 
            col='dataset', 
            kind='kde',
            fill=True,
            height=5,
            aspect=1.2,
            facet_kws={'sharey': False, 'sharex': False} # Un-share axes for better individual plot scaling
        )
        g.fig.suptitle(f'Log-Transformed Distribution for {gene}', fontsize=16, fontweight='bold', y=1.03)
        g.set_axis_labels(f'Log-Transformed {gene} Expression', 'Density')
        g.set_titles("Dataset: {col_name}")
        
        # Remove the temporary column
        combined_df.drop(columns=[log_gene_col], inplace=True)

        filepath = os.path.join(output_dir, f"{gene}_distribution_comparison.png")
        plt.savefig(filepath, bbox_inches='tight')
        plt.close(g.fig)
        print(f"    Saved plot to {filepath}")
        
    print(f"\nAll distribution plots saved to: {output_dir}")


def generate_faceted_plots(datasets, plot_kind, required_genes, output_dir, debug=False):
    """Generate faceted plots (e.g., violin, scatter) for all datasets."""
    print(f"\n{'='*100}")
    print(f"GENERATING {plot_kind.upper()} PLOTS")
    print(f"{'='*100}")
    os.makedirs(output_dir, exist_ok=True)
    
    # Combine data from all datasets for faceted plotting
    combined_df = pd.concat(
        [prepare_data_for_plotting(df, name, debug=debug) for name, df in datasets.items()],
        ignore_index=True
    )
    
    # Check if all required genes for the plot exist in the combined data
    available_genes = [gene for gene in required_genes if gene in combined_df.columns]
    if len(available_genes) != len(required_genes):
        print(f"WARNING: Required genes ({required_genes}) not all present. Skipping {plot_kind} plots.")
        return
        
    print(f"Generating plots for genes: {available_genes}")
    
    if plot_kind == 'violin':
        for gene in available_genes:
            print(f"  - Creating plot for {gene}...")
            # Use catplot for faceted violin plots with split violins
            g = sns.catplot(
                data=combined_df, 
                x='meta_source', 
                y=gene, 
                hue='er_status_str',
                col='dataset', 
                kind='violin',
                split=True, # Create split violins for direct comparison
                inner='quartile',
                height=10, 
                aspect=0.3,
                sharey=False
            )
            g.fig.suptitle(f'Gene Expression for {gene} by Source and ER Status', y=1.03, fontsize=16, fontweight='bold')
            g.set_axis_labels("Source", f"{gene} Expression Level")
            g.set_titles("Dataset: {col_name}")
            g.set_xticklabels(rotation=15, ha='right')
            g.add_legend(title="ER Status")
            
            filepath = os.path.join(output_dir, f"{gene}_split_violin_plot.png")
            plt.savefig(filepath, bbox_inches='tight')
            plt.close('all')
            print(f"    Saved plot to {filepath}")

    elif plot_kind == 'scatter':
        print("  - Creating scatter plot for ESR1 vs. LRRC8D...")
        # Use FacetGrid for customized scatter plots
        g = sns.FacetGrid(
            combined_df, 
            col="dataset", 
            hue="er_status_str", 
            col_wrap=3, 
            height=5,
            hue_order=['Positive', 'Negative']
        )
        g.map(sns.scatterplot, 'ESR1', 'LRRC8D', alpha=0.8, edgecolor='w', s=60)
        g.fig.suptitle('ESR1 vs. LRRC8D Expression by ER Status', y=1.03, fontsize=16, fontweight='bold')
        g.set_axis_labels("ESR1 Expression", "LRRC8D Expression")
        g.set_titles("Dataset: {col_name}")
        g.add_legend(title="ER Status")
        
        filepath = os.path.join(output_dir, "ESR1_vs_LRRC8D_scatter.png")
        plt.savefig(filepath, bbox_inches='tight')
        plt.close('all')
        print(f"    Saved plot to {filepath}")

    print(f"\nAll {plot_kind} plots saved to: {output_dir}")


# --- Main Execution ---
if __name__ == "__main__":
    # Define the directory containing your CSV data files
    DATA_DIRECTORY = '/data/gold/gse_20194_62944/'
    
    # Load all datasets from the directory
    all_datasets = load_datasets(DATA_DIRECTORY)
    
    if all_datasets:
        # --- Run Statistical Analysis ---
        analyzed_genes = {}
        for name, df in all_datasets.items():
            analyzed_genes[name] = analyze_dataframe(df, name)

        # --- Generate Visualizations ---
        # Set debug=True to get detailed output on data processing steps
        debug_mode = True
        
        generate_histograms(all_datasets, output_dir="/outputs/figures/histograms/")
        
        # Generate Violin Plots
        violin_genes = ['ESR1', 'GP2', 'PTK6', 'STX17', 'PPP1R3C', 'SASH1', 'LRRC8D']
        generate_faceted_plots(all_datasets, 'violin', violin_genes, 
                               output_dir="/outputs/figures/violin_plots/", debug=debug_mode)
        
        # Generate Scatter Plots
        scatter_genes = ['ESR1', 'LRRC8D']
        generate_faceted_plots(all_datasets, 'scatter', scatter_genes,
                               output_dir="/outputs/figures/scatter_plots/", debug=debug_mode)

        print("\n\nScript finished. All analyses and visualizations are complete.")
    else:
        print("\nScript finished. No data was processed.")