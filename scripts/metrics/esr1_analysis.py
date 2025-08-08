import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from functools import reduce

# Configure pandas display options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

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
    """Prepare a single dataframe for plotting."""
    df_copy = df.copy()
    df_copy['dataset'] = df_name
    if 'meta_er_status' in df_copy.columns:
        if debug:
            print(f"DEBUG: In prepare_data_for_plotting for {df_name}:")
            print(f"  Unique 'meta_er_status' values BEFORE mapping: {df_copy['meta_er_status'].unique()}")
        
        # Map numeric ER status to string representation
        df_copy['er_status_str'] = df_copy['meta_er_status'].map({1: 'Positive', 0: 'Negative'})
        df_copy.dropna(subset=['er_status_str'], inplace=True)
        df_copy['group'] = df_copy['meta_source'] + '-' + df_copy['er_status_str']
        
        if debug:
            print(f"  Value counts for 'group' AFTER mapping:\n{df_copy['group'].value_counts(dropna=False)}")

    elif debug:
        print(f"DEBUG: 'meta_er_status' column NOT FOUND in {df_name}")

    return df_copy

def generate_histograms(datasets, output_dir="/outputs/figures/histograms/"):
    """Generate log-transformed histograms to compare gene distributions."""
    print(f"\n{'='*100}")
    print("GENERATING HISTOGRAMS")
    print(f"{'='*100}")
    os.makedirs(output_dir, exist_ok=True)

    target_genes = ['ESR1', 'GP2', 'PTK6', 'STX17', 'PPP1R3C', 'SASH1']
    
    # Find common genes across all datasets
    common_genes = list(reduce(set.intersection, [set(df.columns) for df in datasets.values()]))
    genes_to_plot = [gene for gene in target_genes if gene in common_genes]
    
    print(f"Generating histograms for {len(genes_to_plot)} common target genes: {genes_to_plot}")

    # Prepare all datasets for plotting
    prepared_datasets = {name: prepare_data_for_plotting(df, name) for name, df in datasets.items()}
    num_datasets = len(datasets)

    for gene in genes_to_plot:
        print(f"  - Processing {gene}...")
        fig, axes = plt.subplots(num_datasets, 1, figsize=(12, 5 * num_datasets), sharex=True)
        # Ensure axes is always a list for consistent indexing
        if num_datasets == 1:
            axes = [axes]
        fig.suptitle(f'Log-Transformed Distribution for {gene}', fontsize=16, fontweight='bold')

        for ax, (name, df) in zip(axes, prepared_datasets.items()):
            df_copy = df.copy()
            # Apply log transform safely
            df_copy[gene] = np.log1p(df_copy[gene] - df_copy[gene].min())
            
            sns.histplot(data=df_copy, x=gene, hue='meta_source', ax=ax, bins=30, kde=True)
            ax.set_title(f'Dataset: {name}', fontsize=12)
            ax.set_ylabel('Frequency')
        
        axes[-1].set_xlabel(f'Log-Transformed {gene} Expression')
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        filepath = os.path.join(output_dir, f"{gene}_histogram_comparison.png")
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"    Saved histogram to {filepath}")
    print(f"\nAll histograms saved to: {output_dir}")

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
        plot_order = [
            'gse20194-Negative', 'gse20194-Positive',
            'gse62944-Negative', 'gse62944-Positive'
        ]
        palette = {
            'gse20194-Negative': '#3498db', 'gse20194-Positive': '#e74c3c',
            'gse62944-Negative': '#2ecc71', 'gse62944-Positive': '#f1c40f'
        }
        for gene in available_genes:
            print(f"  - Creating plot for {gene}...")
            g = sns.catplot(
                data=combined_df, x='group', y=gene, col='dataset', kind='violin',
                order=plot_order, palette=palette, inner='quartile',
                col_order=list(datasets.keys()), height=10, aspect=(2/len(datasets.keys())),
                sharey=False
            )
            g.fig.suptitle(f'Gene Expression for {gene} by Source and ER Status', y=1.03, fontsize=16, fontweight='bold')
            g.set_axis_labels("Group (Source - ER Status)", f"{gene} Expression Level")
            g.set_titles("Dataset: {col_name}")
            g.set_xticklabels(rotation=45, ha='right')
            filepath = os.path.join(output_dir, f"{gene}_violin_plot.png")
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
            plt.close('all')
            print(f"    Saved plot to {filepath}")

    elif plot_kind == 'scatter':
        print("  - Creating scatter plot for ESR1 vs. LRRC8D...")
        g = sns.FacetGrid(
            combined_df, col="dataset", hue="er_status_str", col_wrap=3, height=5,
            hue_order=['Positive', 'Negative'], palette={'Positive': '#e74c3c', 'Negative': '#3498db'}
        )
        g.map(sns.scatterplot, 'ESR1', 'LRRC8D', alpha=0.7, edgecolor='w', s=50)
        g.fig.suptitle('ESR1 vs. LRRC8D Expression by ER Status', y=1.05, fontsize=16, fontweight='bold')
        g.set_axis_labels("ESR1 Expression", "LRRC8D Expression")
        g.set_titles("Dataset: {col_name}")
        g.add_legend(title="ER Status")
        filepath = os.path.join(output_dir, "ESR1_vs_LRRC8D_scatter.png")
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close('all')
        print(f"    Saved plot to {filepath}")

    print(f"\nAll {plot_kind} plots saved to: {output_dir}")

# --- Main Execution ---
if __name__ == "__main__":
    # Define the directory containing your CSV data files
    DATA_DIRECTORY = '/data/gse_20194_62944/'
    
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