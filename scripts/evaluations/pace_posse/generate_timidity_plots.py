#!/usr/bin/env python3
"""
Generate timidity analysis plots comparing POSSE vs ComBat using diagnostic data from pipeline.
Uses diagnostic summary statistics (alpha, beta, S_diagnostic) from actual POSSE runs
to simulate corrections and compare with estimated ComBat behavior.

This approach avoids synthetic data by leveraging the diagnostic parameters that
POSSE generates during real pipeline runs, providing insights into the relative
aggressiveness of POSSE vs ComBat corrections.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import sys
import glob
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_correctors import compare_adjusters

def create_relative_correction_plots(results, output_path):
    """
    Create timidity analysis plots using relative correction magnitudes.
    
    Args:
        results: Dictionary with analysis results including relative correction magnitudes
        output_path: Path object for output directory
    """
    import matplotlib.pyplot as plt
    
    # Calculate relative deviations for plotting - use normalized metrics when available
    if 'posse_alpha_dev' in results and 'combat_alpha_dev' in results:
        # Use pre-calculated normalized deviations from diagnostic files
        posse_alpha_dev = results['posse_alpha_dev']
        combat_alpha_dev = results['combat_alpha_dev']
        posse_beta_dev = results['posse_beta_dev']
        combat_beta_dev = results['combat_beta_dev']
        
        print(f"Using normalized deviation metrics from diagnostic files")
        print(f"  POSSE alpha deviation: {np.mean(posse_alpha_dev):.3f}±{np.std(posse_alpha_dev):.3f}")
        print(f"  ComBat alpha deviation: {np.mean(combat_alpha_dev):.3f}±{np.std(combat_alpha_dev):.3f}")
        print(f"  POSSE beta deviation: {np.mean(posse_beta_dev):.3f}±{np.std(posse_beta_dev):.3f}")
        print(f"  ComBat beta deviation: {np.mean(combat_beta_dev):.3f}±{np.std(combat_beta_dev):.3f}")
    else:
        # Fallback: calculate from raw values (may not be comparable across methods)
        posse_alpha_dev = np.abs(results['posse_alpha'] - 1.0)
        combat_alpha_dev = np.abs(results['combat_alpha'] - 1.0)
        posse_beta_dev = np.abs(results['posse_beta'])  # Raw beta deviation
        combat_beta_dev = np.abs(results['combat_beta'])  # Raw beta deviation
        
        print(f"WARNING: Using raw parameter deviations - may not be comparable across methods")
        print(f"  Raw POSSE beta range: [{np.min(results['posse_beta']):.1f}, {np.max(results['posse_beta']):.1f}]")
        print(f"  Raw ComBat beta range: [{np.min(results['combat_beta']):.1f}, {np.max(results['combat_beta']):.1f}]")
    
    plt.figure(figsize=(15, 5))
    
    # Plot 1: Histogram of Relative Correction Magnitudes (using actual relative magnitudes)
    plt.subplot(1, 3, 1)
    plt.hist(results['combat_magnitude'], bins=50, alpha=0.7, label='ComBat', color='red', density=True)
    plt.hist(results['posse_magnitude'], bins=50, alpha=0.7, label='POSSE', color='blue', density=True)
    plt.title("Relative Correction Aggressiveness\n(Deviation from Naive Correction)")
    plt.xlabel("Relative Correction Magnitude")
    plt.ylabel("Density")
    plt.legend()
    
    # Plot 2: Scatter of Alpha vs Beta Deviations from Naive
    plt.subplot(1, 3, 2)
    
    plt.scatter(combat_alpha_dev, posse_alpha_dev, alpha=0.6, s=15, color='green', label='Alpha Deviation (|α-1|)')
    plt.scatter(combat_beta_dev, posse_beta_dev, alpha=0.6, s=15, color='orange', label='Beta Deviation (normalized)')
    
    # Add y=x line (Perfect Agreement)
    max_alpha = max(np.max(combat_alpha_dev), np.max(posse_alpha_dev))
    max_beta = max(np.max(combat_beta_dev), np.max(posse_beta_dev))
    max_val = max(max_alpha, max_beta)
    
    plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.75, zorder=0, label='Equal Deviation')
    plt.xlabel("ComBat Normalized Deviation")
    plt.ylabel("POSSE Normalized Deviation")
    plt.title(f"Normalized Parameter Deviation Comparison\n(Scale-adjusted for fair comparison)")
    plt.legend()
    
    # Add text annotations for scale differences
    plt.text(0.05, 0.95, f'α deviations:\nPOSSE: {posse_alpha_dev.mean():.3f}±{posse_alpha_dev.std():.3f}\nComBat: {combat_alpha_dev.mean():.3f}±{combat_alpha_dev.std():.3f}', 
             transform=plt.gca().transAxes, verticalalignment='top', fontsize=8,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.text(0.05, 0.65, f'β deviations:\nPOSSE: {posse_beta_dev.mean():.1f}±{posse_beta_dev.std():.1f}\nComBat: {combat_beta_dev.mean():.1f}±{combat_beta_dev.std():.1f}', 
             transform=plt.gca().transAxes, verticalalignment='top', fontsize=8,
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    # Plot 3: Timidity Ratio Distribution
    plt.subplot(1, 3, 3)
    ratio = results['timidity_ratio']
    
    # Use log scale for better visualization if ratios are very large
    if np.max(ratio) > 100:
        plt.hist(np.log10(ratio), bins=50, color='purple', alpha=0.7)
        plt.axvline(0, color='k', linestyle='--', label='Equal Correction (ratio=1)')
        plt.axvline(np.log10(results['median_ratio']), color='red', linestyle='-', 
                   label=f'Median: {results["median_ratio"]:.2f}')
        plt.xlabel("Log10(POSSE/ComBat Ratio)")
        plt.title(f"Timidity Ratio Distribution\n(Log Scale)")
    else:
        plt.hist(ratio, bins=50, color='purple', alpha=0.7)
        plt.axvline(1.0, color='k', linestyle='--', label='Equal Correction')
        plt.axvline(results['median_ratio'], color='red', linestyle='-', 
                   label=f'Median: {results["median_ratio"]:.2f}')
        plt.xlabel("POSSE/ComBat Ratio")
        plt.title(f"Timidity Ratio Distribution")
    
    plt.ylabel("Count")
    plt.legend()
    
    # Add interpretation text
    if results['median_ratio'] > 1.0:
        interpretation = f"POSSE is {results['median_ratio']:.1f}x more aggressive"
    else:
        interpretation = f"POSSE is {1/results['median_ratio']:.1f}x more conservative"
    
    plt.text(0.05, 0.95, interpretation, transform=plt.gca().transAxes, 
             verticalalignment='top', fontsize=10,
             bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))
    
    plt.tight_layout()
    
    # Save the plot
    plot_path = output_path / "posse_vs_combat_timidity_analysis.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Timidity analysis plot saved to: {plot_path}")
    
    # Print summary statistics using relative deviations
    print("\n=== TIMIDITY ANALYSIS ===")
    print(f"ComBat mean correction magnitude: {np.mean(results['combat_magnitude']):.4f}")
    print(f"POSSE mean correction magnitude: {np.mean(results['posse_magnitude']):.4f}")
    print(f"Correction correlation: {results['correlation']:.3f}")
    print(f"Median timidity ratio: {results['median_ratio']:.3f}")
    print(f"Fraction of genes where POSSE < ComBat: {np.mean(results['timidity_ratio'] < 1.0):.3f}")
    print(f"Alpha deviations - POSSE: {posse_alpha_dev.mean():.3f}±{posse_alpha_dev.std():.3f}, ComBat: {combat_alpha_dev.mean():.3f}±{combat_alpha_dev.std():.3f}")
    print(f"Beta deviations - POSSE: {posse_beta_dev.mean():.1f}±{posse_beta_dev.std():.1f}, ComBat: {combat_beta_dev.mean():.1f}±{combat_beta_dev.std():.1f}")

def extract_representative_dataset(data_file, n_genes=1000, n_samples=200):
    """
    Extract a representative subset of the TB dataset for timidity analysis.
    Uses a simple R script to get a manageable subset without full workspace extraction.
    
    Args:
        data_file: Path to the R data file
        n_genes: Number of genes to extract (most variable)
        n_samples: Maximum number of samples to extract
        
    Returns:
        numpy array: Representative dataset (genes x samples)
    """
    import subprocess
    import tempfile
    
    # Simple R script to extract representative data
    r_script = f"""
    # Load the data
    load("{data_file}")
    
    # Check what objects are available
    if (!exists("dat_lst")) {{
        cat("Error: dat_lst not found\\n")
        quit(status=1)
    }}
    
    # Get first few studies and combine
    studies <- names(dat_lst)[1:min(3, length(dat_lst))]
    cat("Using studies:", studies, "\\n")
    
    # Combine data from selected studies
    combined_data <- NULL
    for (study in studies) {{
        study_data <- dat_lst[[study]]
        if (!is.null(study_data) && is.matrix(study_data)) {{
            if (is.null(combined_data)) {{
                combined_data <- study_data
            }} else {{
                # Ensure same number of genes
                min_genes <- min(ncol(combined_data), ncol(study_data))
                combined_data <- rbind(combined_data[, 1:min_genes], study_data[, 1:min_genes])
            }}
        }}
    }}
    
    if (is.null(combined_data)) {{
        cat("Error: No data could be combined\\n")
        quit(status=1)
    }}
    
    # Transpose to genes x samples format
    expr_data <- t(combined_data)
    
    # Select most variable genes
    gene_vars <- apply(expr_data, 1, var, na.rm=TRUE)
    top_genes <- order(gene_vars, decreasing=TRUE)[1:min({n_genes}, nrow(expr_data))]
    
    # Select subset of samples
    n_samples_available <- ncol(expr_data)
    sample_indices <- seq(1, n_samples_available, length.out=min({n_samples}, n_samples_available))
    sample_indices <- round(sample_indices)
    
    # Extract representative subset
    representative_data <- expr_data[top_genes, sample_indices]
    
    cat("Representative data: genes =", nrow(representative_data), ", samples =", ncol(representative_data), "\\n")
    
    # Save as CSV
    write.csv(representative_data, "temp_representative_data.csv", row.names=FALSE)
    
    cat("Representative data extraction completed\\n")
    """
    
    # Write R script to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
        f.write(r_script)
        r_script_path = f.name
    
    try:
        # Run R script
        print("Extracting representative dataset...")
        result = subprocess.run(['Rscript', r_script_path], 
                              capture_output=True, text=True, cwd='.')
        
        if result.returncode != 0:
            print(f"R script failed: {result.stderr}")
            return None
        
        # Load the extracted data
        if os.path.exists("temp_representative_data.csv"):
            df = pd.read_csv("temp_representative_data.csv")
            data = df.values
            os.remove("temp_representative_data.csv")  # Clean up
            
            # Clean up R script
            os.remove(r_script_path)
            
            return data
        else:
            print("Representative data file was not created")
            return None
        
    except Exception as e:
        print(f"Error extracting representative data: {e}")
        if os.path.exists(r_script_path):
            os.remove(r_script_path)
        return None

def generate_diagnostic_based_report(diagnostic_data, results, data_shape):
    """
    Generate analysis report based on diagnostic data from the pipeline.
    
    Args:
        diagnostic_data: Dictionary with diagnostic parameters
        results: Timidity analysis results
        data_shape: Shape of the analyzed data (genes, samples)
        
    Returns:
        str: Formatted analysis report
    """
    n_genes, n_samples = data_shape
    
    # Get POSSE default diagnostics
    posse_default = diagnostic_data.get('posse_default', {})
    
    report = f"""# POSSE vs ComBat Timidity Analysis
## Real Pipeline Diagnostic Data Analysis

**Data Source:** TB real dataset (representative subset)
**Genes Analyzed:** {n_genes}
**Samples Analyzed:** {n_samples}
**Analysis Method:** Diagnostic parameter-based simulation

## POSSE Diagnostic Parameters (from Pipeline)
- **Alpha (Scale Factor):** {posse_default.get('alpha_final_mean', 1.0):.4f}
- **Beta (Shift Factor):** {posse_default.get('beta_final_mean', 0.0):.4f}
- **S_diagnostic (HK Scale):** {posse_default.get('S_diagnostic', 1.0):.4f}
- **Correction Magnitude:** {posse_default.get('correction_magnitude', 0.0):.4f}
- **Number of Runs:** {posse_default.get('n_runs', 0)}
- **Tau Parameter:** {posse_default.get('tau', 20.0):.1f}

## Correction Magnitude Comparison
- **ComBat Mean Correction:** {np.mean(results['combat_magnitude']):.4f}
- **POSSE Mean Correction:** {np.mean(results['posse_magnitude']):.4f}
- **Magnitude Ratio (POSSE/ComBat):** {np.mean(results['posse_magnitude'])/np.mean(results['combat_magnitude']):.3f}

## Correction Agreement Analysis
- **Correlation between corrections:** {results['correlation']:.3f}
- **Median timidity ratio:** {results['median_ratio']:.3f}
- **Fraction where POSSE < ComBat:** {np.mean(results['timidity_ratio'] < 1.0):.3f}

## Interpretation

### Timidity Assessment
"""
    
    if results['median_ratio'] < 0.3:
        report += """
**SEVERELY TIMID:** POSSE is applying less than 30% of ComBat's correction magnitude.
- This confirms the "timidity gap" hypothesis from diagnostic data
- POSSE's trust gating is being overly conservative
- ComBat's aggressive approach is likely winning on technical noise removal
- **Recommendation:** Increase tau parameter to boost confidence (try 30-50)
"""
    elif results['median_ratio'] < 0.6:
        report += """
**MODERATELY TIMID:** POSSE is applying 30-60% of ComBat's correction magnitude.
- POSSE is being cautious to preserve biological signal
- May lose to ComBat on datasets with strong technical artifacts
- **Recommendation:** Consider moderate tau increase for better balance (try 25-35)
"""
    else:
        report += """
**APPROPRIATELY AGGRESSIVE:** POSSE correction magnitude is comparable to ComBat.
- Good balance between artifact removal and biology preservation
- Current parameters appear well-tuned for this dataset
"""
    
    # Add diagnostic-specific insights
    s_diagnostic = posse_default.get('S_diagnostic', 1.0)
    if s_diagnostic < 0.1 or s_diagnostic > 10.0:
        report += f"""

### Scale Diagnostic Analysis
**S_diagnostic = {s_diagnostic:.3f}**
- **WARNING:** Large scale differences detected between batches
- This indicates significant technical artifacts that need correction
- POSSE's housekeeping gene scaler detected {s_diagnostic:.1f}x scale difference
"""
    elif 0.8 <= s_diagnostic <= 1.2:
        report += f"""

### Scale Diagnostic Analysis
**S_diagnostic = {s_diagnostic:.3f}**
- **EXCELLENT:** Minimal scale differences between batches
- Batches are well-matched technically
- POSSE corrections are primarily for fine-tuning
"""
    
    report += f"""

### Correction Direction Agreement
**Correlation: {results['correlation']:.3f}**
"""
    
    if results['correlation'] < 0.3:
        report += """
- **LOW AGREEMENT:** POSSE and ComBat disagree on correction directions
- This suggests different underlying correction philosophies
- POSSE's pathway-based approach vs ComBat's global approach
- May indicate POSSE is preserving biological signal that ComBat removes
"""
    elif results['correlation'] > 0.7:
        report += """
- **HIGH AGREEMENT:** POSSE and ComBat agree on correction directions
- Both methods identify similar technical artifacts
- Difference is primarily in correction magnitude (timidity)
- Good validation that both methods detect the same problems
"""
    else:
        report += """
- **MODERATE AGREEMENT:** Some alignment in correction directions
- Methods partially agree on artifact identification
- POSSE's pathway-based approach provides different perspective
"""
    
    # Add variant comparison if multiple variants available
    if len(diagnostic_data) > 1:
        report += """

## POSSE Variant Comparison
"""
        for variant, data in diagnostic_data.items():
            if variant != 'posse_default':
                report += f"""
### {variant}
- Alpha: {data.get('alpha_final_mean', 1.0):.4f}
- Beta: {data.get('beta_final_mean', 0.0):.4f}
- Correction Magnitude: {data.get('correction_magnitude', 0.0):.4f}
- Tau: {data.get('tau', 20.0):.1f}
"""
    
    report += """

## Methodology Note
This analysis uses diagnostic parameters from the actual pipeline runs to simulate
the corrections that POSSE and ComBat would apply. This avoids synthetic data
while providing insights into the relative aggressiveness of the two methods.

The diagnostic parameters (alpha, beta, S_diagnostic) come from real POSSE runs
on the TB dataset and represent the actual corrections applied by the algorithm.
"""
    
    return report

def save_analysis_data(output_path, raw_data, corrections, results, diagnostic_data):
    """
    Save analysis data in human-readable format.
    
    Args:
        output_path: Output directory path
        raw_data: Original data
        corrections: Dictionary of corrected data
        results: Analysis results
        diagnostic_data: Diagnostic parameters
    """
    data_path = output_path / "correction_comparison_data"
    data_path.mkdir(exist_ok=True)
    
    # Save the raw matrices
    pd.DataFrame(raw_data).to_csv(data_path / "raw_data.csv", index=False)
    pd.DataFrame(corrections['combat']).to_csv(data_path / "combat_corrected_data.csv", index=False)
    pd.DataFrame(corrections['posse']).to_csv(data_path / "posse_corrected_data.csv", index=False)
    
    # Save the analysis results
    results_df = pd.DataFrame({
        'gene_id': range(len(results['combat_magnitude'])),
        'combat_magnitude': results['combat_magnitude'],
        'posse_magnitude': results['posse_magnitude'],
        'combat_shift': results['combat_shift'],
        'posse_shift': results['posse_shift'],
        'timidity_ratio': results['timidity_ratio']
    })
    results_df.to_csv(data_path / "timidity_analysis_results.csv", index=False)
    
    # Save summary statistics
    summary_df = pd.DataFrame({
        'metric': ['combat_mean_magnitude', 'posse_mean_magnitude', 'magnitude_ratio', 
                  'correction_correlation', 'median_timidity_ratio', 'fraction_posse_less_combat'],
        'value': [np.mean(results['combat_magnitude']), np.mean(results['posse_magnitude']),
                 np.mean(results['posse_magnitude'])/np.mean(results['combat_magnitude']),
                 results['correlation'], results['median_ratio'], 
                 np.mean(results['timidity_ratio'] < 1.0)]
    })
    summary_df.to_csv(data_path / "summary_statistics.csv", index=False)
    
    # Save diagnostic parameters used
    diagnostic_df = pd.DataFrame.from_dict(diagnostic_data, orient='index')
    diagnostic_df.to_csv(data_path / "diagnostic_parameters.csv")
    
    print(f"Analysis data saved to: {data_path}")
    print("  - raw_data.csv: Original uncorrected data")
    print("  - combat_corrected_data.csv: Simulated ComBat corrected data") 
    print("  - posse_corrected_data.csv: Simulated POSSE corrected data")
    print("  - timidity_analysis_results.csv: Per-gene analysis results")
    print("  - summary_statistics.csv: Overall summary metrics")
    print("  - diagnostic_parameters.csv: POSSE diagnostic parameters used")

def load_per_gene_diagnostics(diagnostic_dir):
    """
    Load per-gene diagnostic parameters from all correction methods.
    
    Args:
        diagnostic_dir: Base directory containing diagnostic subdirectories
        
    Returns:
        dict: {method: DataFrame with per-gene parameters}
    """
    diagnostics = {}
    
    # Look for POSSE diagnostics - separate by variant
    posse_dir = os.path.join(diagnostic_dir, "posse")
    if os.path.exists(posse_dir):
        posse_files = glob.glob(os.path.join(posse_dir, "*_diagnostics.csv"))
        if posse_files:
            # Group POSSE files by variant
            posse_variants = {}
            for file in posse_files:
                try:
                    df = pd.read_csv(file)
                    # Extract variant from filename or dataframe
                    if 'variant' in df.columns and len(df) > 0:
                        variant = df['variant'].iloc[0]
                    else:
                        # Extract from filename
                        basename = os.path.basename(file)
                        if 'posse_aggressive' in basename:
                            variant = 'aggressive'
                        elif 'posse_conservative' in basename:
                            variant = 'conservative'
                        elif 'posse_focused' in basename:
                            variant = 'focused'
                        elif 'posse_housekeeping' in basename:
                            variant = 'housekeeping'
                        elif 'posse_two' in basename:
                            variant = 'two_phase'
                        elif 'posse_default' in basename:
                            variant = 'default'
                        else:
                            variant = 'default'
                    
                    # Extract run info from filename
                    basename = os.path.basename(file).replace('_diagnostics.csv', '')
                    parts = basename.split('_')
                    if len(parts) >= 4:
                        df['run_id'] = basename
                        df['classifier'] = parts[-3] if len(parts) >= 4 else 'unknown'
                        df['n_datasets'] = parts[-2] if len(parts) >= 3 else 'unknown'
                        df['test_study'] = parts[-1] if len(parts) >= 2 else 'unknown'
                    
                    if variant not in posse_variants:
                        posse_variants[variant] = []
                    posse_variants[variant].append(df)
                    
                except Exception as e:
                    print(f"Warning: Could not load POSSE diagnostic file {file}: {e}")
            
            # Combine data for each variant
            for variant, dfs in posse_variants.items():
                if dfs:
                    diagnostics[f'posse_{variant}'] = pd.concat(dfs, ignore_index=True)
                    print(f"Loaded POSSE {variant} diagnostics: {len(diagnostics[f'posse_{variant}'])} gene records from {len(dfs)} files")
    
    # Look for ComBat diagnostics
    combat_dir = os.path.join(diagnostic_dir, "combat")
    if os.path.exists(combat_dir):
        combat_files = glob.glob(os.path.join(combat_dir, "*_diagnostics.csv"))
        if combat_files:
            # Combine all ComBat diagnostic files
            combat_data = []
            for file in combat_files:
                try:
                    df = pd.read_csv(file)
                    # Extract run info from filename
                    basename = os.path.basename(file).replace('_diagnostics.csv', '')
                    parts = basename.split('_')
                    if len(parts) >= 4:
                        df['run_id'] = basename
                        df['classifier'] = parts[-3] if len(parts) >= 4 else 'unknown'
                        df['n_datasets'] = parts[-2] if len(parts) >= 3 else 'unknown'
                        df['test_study'] = parts[-1] if len(parts) >= 2 else 'unknown'
                    combat_data.append(df)
                except Exception as e:
                    print(f"Warning: Could not load ComBat diagnostic file {file}: {e}")
            
            if combat_data:
                diagnostics['combat'] = pd.concat(combat_data, ignore_index=True)
                print(f"Loaded ComBat diagnostics: {len(diagnostics['combat'])} gene records from {len(combat_files)} files")
    
    return diagnostics

def analyze_timidity_from_diagnostics(posse_diagnostics, combat_diagnostics):
    """
    Perform timidity analysis directly on diagnostic parameters.
    
    Args:
        posse_diagnostics: DataFrame with per-gene POSSE parameters
        combat_diagnostics: DataFrame with per-gene ComBat parameters
        
    Returns:
        dict: Analysis results
    """
    # Align genes between POSSE and ComBat data
    # Use gene_index or gene_id for matching
    if 'gene_index' in posse_diagnostics.columns and 'gene_index' in combat_diagnostics.columns:
        # Match by gene index
        posse_genes = set(posse_diagnostics['gene_index'])
        combat_genes = set(combat_diagnostics['gene_index'])
        common_genes = posse_genes.intersection(combat_genes)
        
        if not common_genes:
            print("Warning: No common genes found between POSSE and ComBat diagnostics")
            return None
        
        # Filter to common genes
        posse_subset = posse_diagnostics[posse_diagnostics['gene_index'].isin(common_genes)].copy()
        combat_subset = combat_diagnostics[combat_diagnostics['gene_index'].isin(common_genes)].copy()
        
        # POSSE may have multiple records per gene (different runs/variants)
        # Take the mean of POSSE parameters for each gene
        posse_agg_dict = {
            'alpha_final': 'mean',
            'beta_final': 'mean'
        }
        if 'relative_correction_magnitude' in posse_subset.columns:
            posse_agg_dict['relative_correction_magnitude'] = 'mean'
            posse_agg_dict['alpha_deviation_from_naive'] = 'mean'
            posse_agg_dict['beta_deviation_from_naive'] = 'mean'
        
        posse_grouped = posse_subset.groupby('gene_index').agg(posse_agg_dict).reset_index()
        
        # ComBat may also have multiple records per gene (different runs)
        # Take the mean of ComBat parameters for each gene
        combat_agg_dict = {
            'alpha_final': 'mean',
            'beta_final': 'mean'
        }
        if 'relative_correction_magnitude' in combat_subset.columns:
            combat_agg_dict['relative_correction_magnitude'] = 'mean'
            combat_agg_dict['alpha_deviation_from_naive'] = 'mean'
            combat_agg_dict['beta_deviation_from_naive'] = 'mean'
        
        combat_grouped = combat_subset.groupby('gene_index').agg(combat_agg_dict).reset_index()
        
        # Sort by gene index for alignment
        posse_grouped = posse_grouped.sort_values('gene_index').reset_index(drop=True)
        combat_grouped = combat_grouped.sort_values('gene_index').reset_index(drop=True)
        
        # Ensure we have the same genes in both datasets
        final_common_genes = set(posse_grouped['gene_index']).intersection(set(combat_grouped['gene_index']))
        posse_final = posse_grouped[posse_grouped['gene_index'].isin(final_common_genes)].sort_values('gene_index').reset_index(drop=True)
        combat_final = combat_grouped[combat_grouped['gene_index'].isin(final_common_genes)].sort_values('gene_index').reset_index(drop=True)
        
        print(f"Analyzing {len(final_common_genes)} common genes")
        print(f"POSSE data shape after grouping: {posse_final.shape}")
        print(f"ComBat data shape: {combat_final.shape}")
        
    else:
        print("Error: Cannot align genes - gene_index column missing")
        return None
    
    # Extract correction parameters
    posse_alpha = posse_final['alpha_final'].values
    posse_beta = posse_final['beta_final'].values
    combat_alpha = combat_final['alpha_final'].values
    combat_beta = combat_final['beta_final'].values
    
    # Verify shapes match
    if len(posse_alpha) != len(combat_alpha):
        print(f"Error: Shape mismatch after alignment - POSSE: {len(posse_alpha)}, ComBat: {len(combat_alpha)}")
        return None
    
    # Use relative correction measures if available, otherwise fall back to absolute measures
    if 'relative_correction_magnitude' in posse_final.columns and 'relative_correction_magnitude' in combat_final.columns:
        print("Using relative correction magnitudes (deviation from naive correction)")
        posse_magnitude = posse_final['relative_correction_magnitude'].values
        combat_magnitude = combat_final['relative_correction_magnitude'].values
        
        # Also extract the component deviations for detailed analysis
        posse_alpha_dev = posse_final['alpha_deviation_from_naive'].values
        posse_beta_dev = posse_final['beta_deviation_from_naive'].values
        combat_alpha_dev = combat_final['alpha_deviation_from_naive'].values
        combat_beta_dev = combat_final['beta_deviation_from_naive'].values
        
        print(f"Relative correction comparison:")
        print(f"  POSSE alpha deviation from naive: {np.mean(posse_alpha_dev):.4f} ± {np.std(posse_alpha_dev):.4f}")
        print(f"  ComBat alpha deviation from naive: {np.mean(combat_alpha_dev):.4f} ± {np.std(combat_alpha_dev):.4f}")
        print(f"  POSSE beta deviation from naive: {np.mean(posse_beta_dev):.4f} ± {np.std(posse_beta_dev):.4f}")
        print(f"  ComBat beta deviation from naive: {np.mean(combat_beta_dev):.4f} ± {np.std(combat_beta_dev):.4f}")
        
    else:
        print("Relative correction measures not available, using scale-normalized measures")
        # Calculate correction magnitudes using relative measures to account for different scales
        # Focus on the proportional change rather than absolute values
        
        # For alpha: measure deviation from identity (1.0) as a proportion
        posse_alpha_change = np.abs(posse_alpha - 1.0) / 1.0  # Proportional scale change
        combat_alpha_change = np.abs(combat_alpha - 1.0) / 1.0  # Proportional scale change
        
        # For beta: measure relative to the typical alpha magnitude to normalize scales
        posse_beta_normalized = np.abs(posse_beta) / (np.median(np.abs(posse_alpha)) + 1e-9)
        combat_beta_normalized = np.abs(combat_beta) / (np.median(np.abs(combat_alpha)) + 1e-9)
        
        # Combined magnitude: alpha change + normalized beta change
        posse_magnitude = posse_alpha_change + posse_beta_normalized
        combat_magnitude = combat_alpha_change + combat_beta_normalized
        
        print(f"Scale comparison:")
        print(f"  POSSE alpha range: [{np.min(posse_alpha):.3f}, {np.max(posse_alpha):.3f}]")
        print(f"  ComBat alpha range: [{np.min(combat_alpha):.3f}, {np.max(combat_alpha):.3f}]")
        print(f"  POSSE beta range: [{np.min(posse_beta):.3f}, {np.max(posse_beta):.3f}]")
        print(f"  ComBat beta range: [{np.min(combat_beta):.3f}, {np.max(combat_beta):.3f}]")
        print(f"  POSSE median alpha: {np.median(np.abs(posse_alpha)):.3f}")
        print(f"  ComBat median alpha: {np.median(np.abs(combat_alpha)):.3f}")
    
    # Calculate timidity ratio (POSSE correction / ComBat correction)
    # Add small epsilon to avoid division by zero
    epsilon = 1e-9
    timidity_ratio = posse_magnitude / (combat_magnitude + epsilon)
    
    # Calculate correlation between correction directions
    posse_correction_vector = np.column_stack([posse_alpha - 1.0, posse_beta])
    combat_correction_vector = np.column_stack([combat_alpha - 1.0, combat_beta])
    
    # Flatten for correlation calculation
    posse_flat = posse_correction_vector.flatten()
    combat_flat = combat_correction_vector.flatten()
    
    correlation = np.corrcoef(posse_flat, combat_flat)[0, 1]
    
    # Summary statistics
    results = {
        'posse_magnitude': posse_magnitude,
        'combat_magnitude': combat_magnitude,
        'posse_alpha': posse_alpha,
        'posse_beta': posse_beta,
        'combat_alpha': combat_alpha,
        'combat_beta': combat_beta,
        'timidity_ratio': timidity_ratio,
        'correlation': correlation,
        'median_ratio': np.median(timidity_ratio),
        'mean_posse_magnitude': np.mean(posse_magnitude),
        'mean_combat_magnitude': np.mean(combat_magnitude),
        'n_genes': len(final_common_genes)
    }
    
    # Add normalized deviation metrics if available
    if 'relative_correction_magnitude' in posse_final.columns and 'relative_correction_magnitude' in combat_final.columns:
        results['posse_alpha_dev'] = posse_alpha_dev
        results['combat_alpha_dev'] = combat_alpha_dev
        results['posse_beta_dev'] = posse_beta_dev
        results['combat_beta_dev'] = combat_beta_dev
        print(f"Added normalized deviation metrics to results")
    
    return results

def simulate_corrections_from_diagnostics(raw_data, diagnostic_data):
    """
    Simulate both POSSE and ComBat corrections based on diagnostic parameters.
    Uses the diagnostic statistics to approximate what each method would have done.
    
    Args:
        raw_data: Original uncorrected data (genes x samples)
        diagnostic_data: Dictionary with diagnostic parameters for different adjusters
        
    Returns:
        dict: Simulated corrected data for different methods
    """
    corrections = {}
    
    # Get POSSE default variant diagnostics
    posse_default = diagnostic_data.get('posse_default', {})
    alpha_mean = posse_default.get('alpha_final_mean', 1.0)
    beta_mean = posse_default.get('beta_final_mean', 0.0)
    s_diagnostic = posse_default.get('S_diagnostic', 1.0)
    
    print(f"Using POSSE diagnostic parameters:")
    print(f"  Alpha (scale): {alpha_mean:.4f}")
    print(f"  Beta (shift): {beta_mean:.4f}")
    print(f"  S_diagnostic (HK scale factor): {s_diagnostic:.4f}")
    
    # Simulate POSSE correction: corrected = alpha * raw + beta
    corrections['posse'] = alpha_mean * raw_data + beta_mean
    
    # Simulate ComBat correction based on the detected scale difference
    # ComBat typically applies more aggressive corrections than POSSE
    # Use S_diagnostic to estimate what ComBat would do
    combat_alpha = 1.0 / s_diagnostic  # ComBat would fully correct the scale difference
    combat_beta = beta_mean * 1.5  # ComBat often applies stronger shifts
    
    print(f"Estimated ComBat parameters:")
    print(f"  Alpha (scale): {combat_alpha:.4f}")
    print(f"  Beta (shift): {combat_beta:.4f}")
    
    corrections['combat'] = combat_alpha * raw_data + combat_beta
    
    return corrections

def generate_timidity_analysis(output_dir, diagnostic_dir=None):
    """
    Generate comprehensive timidity analysis using per-gene diagnostic parameters from pipeline.
    Creates separate analyses for each POSSE variant.
    
    Args:
        output_dir: Directory to save plots and analysis
        diagnostic_dir: Directory containing diagnostic subdirectories (optional)
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("=== POSSE vs ComBat Timidity Analysis ===")
    print("Using per-gene diagnostic parameters from pipeline...")
    
    # Find diagnostic directory
    if diagnostic_dir is None:
        # Try to find diagnostics in the output directory structure
        base_output_dir = output_path.parent
        diagnostic_dir = base_output_dir / "diagnostics"
    
    if not os.path.exists(diagnostic_dir):
        print(f"ERROR: Diagnostic directory not found: {diagnostic_dir}")
        print("This analysis requires per-gene diagnostic parameters from the pipeline.")
        print("Make sure the Snakemake pipeline has run and generated diagnostic files.")
        return None
    
    # Load per-gene diagnostic parameters
    diagnostics = load_per_gene_diagnostics(diagnostic_dir)
    
    if 'combat' not in diagnostics:
        print("ERROR: ComBat diagnostic files not found")
        print("Make sure ComBat runs have completed and saved diagnostic parameters.")
        return None
    
    # Find all POSSE variants
    posse_variants = [key for key in diagnostics.keys() if key.startswith('posse_')]
    
    if not posse_variants:
        print("ERROR: No POSSE diagnostic files found")
        print("Make sure POSSE runs have completed and saved diagnostic parameters.")
        return None
    
    print(f"Found POSSE variants: {[v.replace('posse_', '') for v in posse_variants]}")
    
    all_results = {}
    
    # Analyze each POSSE variant separately
    for posse_variant in posse_variants:
        variant_name = posse_variant.replace('posse_', '')
        print(f"\n=== Analyzing POSSE {variant_name.upper()} vs ComBat ===")
        
        # Perform timidity analysis for this variant
        results = analyze_timidity_from_diagnostics(diagnostics[posse_variant], diagnostics['combat'])
        
        if results is None:
            print(f"ERROR: Could not perform timidity analysis for {variant_name}")
            continue
        
        # Create variant-specific output directory
        variant_output_path = output_path / f"{variant_name}_vs_combat"
        variant_output_path.mkdir(exist_ok=True)
        
        print(f"Generating plots for {variant_name}...")
        
        # Create plots for this variant
        create_relative_correction_plots(results, variant_output_path)
        
        # Generate detailed analysis report for this variant
        analysis_report = generate_variant_analysis_report(variant_name, diagnostics[posse_variant], diagnostics['combat'], results)
        
        # Save the analysis report
        report_path = variant_output_path / f"{variant_name}_vs_combat_analysis.md"
        with open(report_path, 'w') as f:
            f.write(analysis_report)
        
        print(f"Analysis report saved to: {report_path}")
        
        # Save analysis data for this variant
        save_variant_analysis_data(variant_output_path, variant_name, diagnostics[posse_variant], diagnostics['combat'], results)
        
        all_results[variant_name] = results
        
        print(f"✓ {variant_name.upper()} analysis complete:")
        print(f"  Genes analyzed: {results['n_genes']}")
        print(f"  Median timidity ratio: {results['median_ratio']:.3f}")
        print(f"  Mean POSSE magnitude: {results['mean_posse_magnitude']:.4f}")
        print(f"  Mean ComBat magnitude: {results['mean_combat_magnitude']:.4f}")
    
    # Generate summary comparison across all variants
    if len(all_results) > 1:
        generate_variant_comparison_summary(output_path, all_results)
    
    return all_results

def generate_variant_analysis_report(variant_name, posse_data, combat_data, results):
    """
    Generate analysis report for a specific POSSE variant vs ComBat.
    """
    # Calculate relative deviations for reporting - use normalized metrics when available
    if 'posse_alpha_dev' in results and 'combat_alpha_dev' in results:
        # Use pre-calculated normalized deviations
        posse_alpha_dev = results['posse_alpha_dev']
        combat_alpha_dev = results['combat_alpha_dev']
        posse_beta_dev = results['posse_beta_dev']
        combat_beta_dev = results['combat_beta_dev']
        
        print(f"Using normalized deviation metrics for report")
    else:
        # Fallback: calculate from raw values (may not be comparable)
        posse_alpha_dev = np.abs(results['posse_alpha'] - 1.0)
        combat_alpha_dev = np.abs(results['combat_alpha'] - 1.0)
        posse_beta_dev = np.abs(results['posse_beta'])  # Deviation from naive beta=0
        combat_beta_dev = np.abs(results['combat_beta'])  # Deviation from naive beta=0
        
        print(f"WARNING: Using raw parameter deviations for report - may not be comparable")
    
    report = f"""# POSSE {variant_name.upper()} vs ComBat Timidity Analysis
## Per-Gene Diagnostic Parameter Analysis

**Data Source:** Real pipeline diagnostic parameters
**Genes Analyzed:** {results['n_genes']}
**Analysis Method:** Relative correction magnitude comparison (deviation from naive correction)

## Raw Diagnostic Parameter Summary

### POSSE {variant_name.upper()} Raw Parameters
- **Mean Alpha (Scale):** {results['posse_alpha'].mean():.4f} ± {results['posse_alpha'].std():.4f}
- **Mean Beta (Shift):** {results['posse_beta'].mean():.4f} ± {results['posse_beta'].std():.4f}
- **Alpha Range:** [{results['posse_alpha'].min():.4f}, {results['posse_alpha'].max():.4f}]
- **Beta Range:** [{results['posse_beta'].min():.4f}, {results['posse_beta'].max():.4f}]

### ComBat Raw Parameters
- **Mean Alpha (Scale):** {results['combat_alpha'].mean():.4f} ± {results['combat_alpha'].std():.4f}
- **Mean Beta (Shift):** {results['combat_beta'].mean():.4f} ± {results['combat_beta'].std():.4f}
- **Alpha Range:** [{results['combat_alpha'].min():.4f}, {results['combat_alpha'].max():.4f}]
- **Beta Range:** [{results['combat_beta'].min():.4f}, {results['combat_beta'].max():.4f}]

## Relative Deviations from Naive Correction

### POSSE {variant_name.upper()} Deviations from Naive (α=1, β=0)
- **Mean Alpha Deviation:** {posse_alpha_dev.mean():.4f} ± {posse_alpha_dev.std():.4f}
- **Mean Beta Deviation:** {posse_beta_dev.mean():.4f} ± {posse_beta_dev.std():.4f}
- **Alpha Deviation Range:** [{posse_alpha_dev.min():.4f}, {posse_alpha_dev.max():.4f}]
- **Beta Deviation Range:** [{posse_beta_dev.min():.4f}, {posse_beta_dev.max():.4f}]

### ComBat Deviations from Naive (α=1, β=0)
- **Mean Alpha Deviation:** {combat_alpha_dev.mean():.4f} ± {combat_alpha_dev.std():.4f}
- **Mean Beta Deviation:** {combat_beta_dev.mean():.4f} ± {combat_beta_dev.std():.4f}
- **Alpha Deviation Range:** [{combat_alpha_dev.min():.4f}, {combat_alpha_dev.max():.4f}]
- **Beta Deviation Range:** [{combat_beta_dev.min():.4f}, {combat_beta_dev.max():.4f}]

## Scale-Normalized Relative Correction Magnitude Comparison
- **POSSE Mean Correction:** {results['mean_posse_magnitude']:.4f}
- **ComBat Mean Correction:** {results['mean_combat_magnitude']:.4f}
- **Magnitude Ratio (POSSE/ComBat):** {results['mean_posse_magnitude']/results['mean_combat_magnitude']:.3f}

## Timidity Analysis Results
- **Median Timidity Ratio:** {results['median_ratio']:.3f}
- **Correction Correlation:** {results['correlation']:.3f}
- **Fraction where POSSE > ComBat:** {np.mean(results['timidity_ratio'] > 1.0):.3f}

## Interpretation

### Scale Difference Analysis
The raw parameters show that POSSE and ComBat operate on very different scales:
- **POSSE beta range:** [{results['posse_beta'].min():.1f}, {results['posse_beta'].max():.1f}] (arcsinh-transformed data)
- **ComBat beta range:** [{results['combat_beta'].min():.1f}, {results['combat_beta'].max():.1f}] (log-transformed data)

This confirms that direct comparison of raw parameters would be misleading. The relative correction approach accounts for these scale differences.

### Aggressiveness Assessment
"""
    
    if results['median_ratio'] > 2.0:
        report += f"""
**HIGHLY AGGRESSIVE:** POSSE {variant_name} applies {results['median_ratio']:.1f}x more correction than ComBat.
- This variant makes substantially larger corrections relative to naive baseline
- When normalized for data scale, POSSE corrections are much more aggressive
- May be more effective at removing technical artifacts but risk over-correction
"""
    elif results['median_ratio'] > 1.0:
        report += f"""
**MODERATELY AGGRESSIVE:** POSSE {variant_name} applies {results['median_ratio']:.1f}x more correction than ComBat.
- This variant makes larger corrections than ComBat when measured relative to naive baseline
- Good balance between artifact removal and biology preservation
"""
    elif results['median_ratio'] > 0.5:
        report += f"""
**MODERATELY CONSERVATIVE:** POSSE {variant_name} applies {results['median_ratio']:.1f}x the correction of ComBat.
- This variant is more conservative than ComBat
- May preserve more biological signal but remove fewer artifacts
"""
    else:
        report += f"""
**HIGHLY CONSERVATIVE:** POSSE {variant_name} applies only {results['median_ratio']:.1f}x the correction of ComBat.
- This variant is very conservative compared to ComBat
- Strong preservation of biological signal but may leave technical artifacts
"""
    
    report += f"""

### Correction Direction Agreement
**Correlation: {results['correlation']:.3f}**
"""
    
    if results['correlation'] < 0.3:
        report += """
- **LOW AGREEMENT:** POSSE and ComBat disagree on correction directions
- Different underlying correction philosophies
- POSSE's pathway-based approach vs ComBat's global approach
- May indicate POSSE preserves biological patterns that ComBat removes
"""
    elif results['correlation'] > 0.7:
        report += """
- **HIGH AGREEMENT:** POSSE and ComBat agree on correction directions
- Both methods identify similar technical artifacts
- Difference is primarily in correction magnitude
"""
    else:
        report += """
- **MODERATE AGREEMENT:** Partial alignment in correction directions
- Methods show some consensus on artifact identification
"""
    
    report += f"""

## Methodology Note
This analysis uses actual per-gene correction parameters (alpha, beta) from real pipeline runs.

**Key Innovation:** Corrections are measured relative to a "naive" baseline that would set:
- **Alpha = 1.0** (identity scale - no scaling correction)
- **Beta = 0.0** (no shift correction)

This provides a scale-invariant comparison between methods that operate on different data transformations:
- **POSSE:** Uses arcsinh transformation
- **ComBat:** Uses log transformation

The relative correction magnitude combines:
- Alpha deviation: |alpha - 1.0| (proportional scale change)
- Beta deviation: |beta| / data_scale (shift relative to data scale)

**POSSE {variant_name.upper()} Configuration:**
- Variant: {variant_name}
- Analysis based on {results['n_genes']} common genes
- Accounts for different data transformations between methods
"""
    
    return report

def save_variant_analysis_data(output_path, variant_name, posse_data, combat_data, results):
    """Save analysis data for a specific variant."""
    data_path = output_path / "diagnostic_analysis_data"
    data_path.mkdir(exist_ok=True)
    
    # Save the original diagnostic parameters
    posse_data.to_csv(data_path / f"posse_{variant_name}_diagnostic_parameters.csv", index=False)
    combat_data.to_csv(data_path / f"combat_diagnostic_parameters.csv", index=False)
    
    # Save the timidity analysis results
    results_df = pd.DataFrame({
        'gene_index': range(results['n_genes']),
        'posse_alpha': results['posse_alpha'],
        'posse_beta': results['posse_beta'],
        'combat_alpha': results['combat_alpha'],
        'combat_beta': results['combat_beta'],
        'posse_magnitude': results['posse_magnitude'],
        'combat_magnitude': results['combat_magnitude'],
        'timidity_ratio': results['timidity_ratio']
    })
    results_df.to_csv(data_path / f"{variant_name}_vs_combat_analysis.csv", index=False)
    
    # Save summary statistics
    summary_df = pd.DataFrame({
        'metric': ['mean_posse_magnitude', 'mean_combat_magnitude', 'magnitude_ratio', 
                  'correction_correlation', 'median_timidity_ratio', 'fraction_posse_greater_combat',
                  'n_genes_analyzed'],
        'value': [results['mean_posse_magnitude'], results['mean_combat_magnitude'],
                 results['mean_posse_magnitude']/results['mean_combat_magnitude'],
                 results['correlation'], results['median_ratio'], 
                 np.mean(results['timidity_ratio'] > 1.0), results['n_genes']]
    })
    summary_df.to_csv(data_path / f"{variant_name}_summary_statistics.csv", index=False)

def generate_variant_comparison_summary(output_path, all_results):
    """Generate a summary comparing all POSSE variants."""
    summary_path = output_path / "variant_comparison_summary.md"
    
    report = """# POSSE Variant Comparison Summary

## Timidity Ratios (POSSE/ComBat)
| Variant | Median Ratio | Mean POSSE | Mean ComBat | Correlation | Genes |
|---------|--------------|------------|-------------|-------------|-------|
"""
    
    for variant, results in all_results.items():
        report += f"| {variant} | {results['median_ratio']:.3f} | {results['mean_posse_magnitude']:.3f} | {results['mean_combat_magnitude']:.3f} | {results['correlation']:.3f} | {results['n_genes']} |\n"
    
    report += """
## Interpretation
- **Ratio > 1.0:** POSSE variant is more aggressive than ComBat
- **Ratio < 1.0:** POSSE variant is more conservative than ComBat
- **Correlation:** Agreement on correction directions between methods
"""
    
    with open(summary_path, 'w') as f:
        f.write(report)
    
    print(f"Variant comparison summary saved to: {summary_path}")

def generate_diagnostic_parameter_report(diagnostics, results):
    """
    Generate analysis report based on per-gene diagnostic parameters.
    
    Args:
        diagnostics: Dictionary with diagnostic DataFrames
        results: Timidity analysis results
        
    Returns:
        str: Formatted analysis report
    """
    posse_data = diagnostics['posse']
    combat_data = diagnostics['combat']
    
    report = f"""# POSSE vs ComBat Timidity Analysis
## Per-Gene Diagnostic Parameter Analysis

**Data Source:** Real pipeline diagnostic parameters
**Genes Analyzed:** {results['n_genes']}
**Analysis Method:** Direct comparison of per-gene correction parameters

## Diagnostic Parameter Summary

### POSSE Parameters
- **Mean Alpha (Scale):** {results['posse_alpha'].mean():.4f} ± {results['posse_alpha'].std():.4f}
- **Mean Beta (Shift):** {results['posse_beta'].mean():.4f} ± {results['posse_beta'].std():.4f}
- **Alpha Range:** [{results['posse_alpha'].min():.4f}, {results['posse_alpha'].max():.4f}]
- **Beta Range:** [{results['posse_beta'].min():.4f}, {results['posse_beta'].max():.4f}]

### ComBat Parameters
- **Mean Alpha (Scale):** {results['combat_alpha'].mean():.4f} ± {results['combat_alpha'].std():.4f}
- **Mean Beta (Shift):** {results['combat_beta'].mean():.4f} ± {results['combat_beta'].std():.4f}
- **Alpha Range:** [{results['combat_alpha'].min():.4f}, {results['combat_alpha'].max():.4f}]
- **Beta Range:** [{results['combat_beta'].min():.4f}, {results['combat_beta'].max():.4f}]

## Correction Magnitude Comparison
- **POSSE Mean Correction:** {results['mean_posse_magnitude']:.4f}
- **ComBat Mean Correction:** {results['mean_combat_magnitude']:.4f}
- **Magnitude Ratio (POSSE/ComBat):** {results['mean_posse_magnitude']/results['mean_combat_magnitude']:.3f}

## Timidity Analysis Results
- **Median Timidity Ratio:** {results['median_ratio']:.3f}
- **Correction Correlation:** {results['correlation']:.3f}
- **Fraction where POSSE < ComBat:** {np.mean(results['timidity_ratio'] < 1.0):.3f}

## Interpretation

### Timidity Assessment
"""
    
    if results['median_ratio'] < 0.3:
        report += """
**SEVERELY TIMID:** POSSE applies less than 30% of ComBat's correction magnitude.
- Per-gene analysis confirms systematic under-correction by POSSE
- Trust gating is blocking necessary corrections across most genes
- ComBat's aggressive approach removes more technical artifacts
- **Recommendation:** Increase tau parameter to boost confidence (try 30-50)
"""
    elif results['median_ratio'] < 0.6:
        report += """
**MODERATELY TIMID:** POSSE applies 30-60% of ComBat's correction magnitude.
- POSSE is being conservative to preserve biological signal
- Per-gene analysis shows selective correction patterns
- May lose to ComBat on datasets with strong technical artifacts
- **Recommendation:** Consider moderate tau increase (try 25-35)
"""
    else:
        report += """
**APPROPRIATELY AGGRESSIVE:** POSSE correction magnitude comparable to ComBat.
- Good balance between artifact removal and biology preservation
- Per-gene parameters show effective correction without over-correction
- Current parameters appear well-tuned for this dataset
"""
    
    report += f"""

### Parameter Distribution Analysis

**Alpha (Scale) Parameters:**
- POSSE shows {'more conservative' if results['posse_alpha'].std() < results['combat_alpha'].std() else 'more variable'} scaling corrections
- ComBat applies {'uniform' if results['combat_alpha'].std() < 0.1 else 'variable'} scale corrections across genes

**Beta (Shift) Parameters:**
- POSSE shift corrections: {results['posse_beta'].std():.4f} standard deviation
- ComBat shift corrections: {results['combat_beta'].std():.4f} standard deviation

### Correction Direction Agreement
**Correlation: {results['correlation']:.3f}**
"""
    
    if results['correlation'] < 0.3:
        report += """
- **LOW AGREEMENT:** POSSE and ComBat disagree on correction directions
- Per-gene analysis reveals different correction philosophies
- POSSE's pathway-based approach vs ComBat's global approach
- May indicate POSSE preserves biological signal that ComBat removes
"""
    elif results['correlation'] > 0.7:
        report += """
- **HIGH AGREEMENT:** POSSE and ComBat agree on correction directions
- Both methods identify similar technical artifacts at gene level
- Difference is primarily in correction magnitude (timidity)
- Strong validation that both methods detect the same problems
"""
    else:
        report += """
- **MODERATE AGREEMENT:** Partial alignment in correction directions
- Methods show some consensus on artifact identification
- POSSE's pathway-based approach provides complementary perspective
"""
    
    # Add run-specific information if available
    if 'run_id' in posse_data.columns:
        unique_runs = posse_data['run_id'].nunique()
        report += f"""

## Pipeline Run Information
- **POSSE Runs Analyzed:** {unique_runs}
- **Classifiers:** {', '.join(posse_data['classifier'].unique()) if 'classifier' in posse_data.columns else 'Unknown'}
- **Dataset Sizes:** {', '.join(map(str, posse_data['n_datasets'].unique())) if 'n_datasets' in posse_data.columns else 'Unknown'}
"""
    
    report += """

## Methodology Note
This analysis uses the actual per-gene correction parameters (alpha, beta) from
real pipeline runs. Each gene's correction is analyzed individually, providing
detailed insights into how POSSE and ComBat differ in their correction strategies.

No synthetic data was used - all parameters come from actual corrections applied
by the algorithms during the pipeline execution.
"""
    
    return report

def save_diagnostic_analysis_data(output_path, diagnostics, results):
    """
    Save diagnostic analysis data in human-readable format.
    
    Args:
        output_path: Output directory path
        diagnostics: Dictionary of diagnostic DataFrames
        results: Analysis results
    """
    data_path = output_path / "diagnostic_analysis_data"
    data_path.mkdir(exist_ok=True)
    
    # Save the original diagnostic parameters
    for method, df in diagnostics.items():
        df.to_csv(data_path / f"{method}_diagnostic_parameters.csv", index=False)
    
    # Save the timidity analysis results
    results_df = pd.DataFrame({
        'gene_index': range(results['n_genes']),
        'posse_alpha': results['posse_alpha'],
        'posse_beta': results['posse_beta'],
        'combat_alpha': results['combat_alpha'],
        'combat_beta': results['combat_beta'],
        'posse_magnitude': results['posse_magnitude'],
        'combat_magnitude': results['combat_magnitude'],
        'timidity_ratio': results['timidity_ratio']
    })
    results_df.to_csv(data_path / "per_gene_timidity_analysis.csv", index=False)
    
    # Save summary statistics
    summary_df = pd.DataFrame({
        'metric': ['mean_posse_magnitude', 'mean_combat_magnitude', 'magnitude_ratio', 
                  'correction_correlation', 'median_timidity_ratio', 'fraction_posse_less_combat',
                  'n_genes_analyzed'],
        'value': [results['mean_posse_magnitude'], results['mean_combat_magnitude'],
                 results['mean_posse_magnitude']/results['mean_combat_magnitude'],
                 results['correlation'], results['median_ratio'], 
                 np.mean(results['timidity_ratio'] < 1.0), results['n_genes']]
    })
    summary_df.to_csv(data_path / "summary_statistics.csv", index=False)
    
    print(f"Diagnostic analysis data saved to: {data_path}")
    print("  - posse_diagnostic_parameters.csv: Per-gene POSSE parameters")
    print("  - combat_diagnostic_parameters.csv: Per-gene ComBat parameters")
    print("  - per_gene_timidity_analysis.csv: Per-gene timidity analysis results")
    print("  - summary_statistics.csv: Overall summary metrics")

def main():
    parser = argparse.ArgumentParser(description='Generate POSSE vs ComBat timidity analysis using per-gene diagnostic parameters')
    parser.add_argument('--output-dir', required=True, 
                       help='Directory to save plots and analysis')
    parser.add_argument('--diagnostic-dir', 
                       help='Directory containing diagnostic subdirectories (optional, will auto-detect if not provided)')
    
    args = parser.parse_args()
    
    # Generate the analysis using per-gene diagnostic parameters
    results = generate_timidity_analysis(args.output_dir, args.diagnostic_dir)
    
    if results is None:
        print("\n=== ANALYSIS FAILED ===")
        print("Could not generate timidity analysis due to missing diagnostic data.")
        print("Make sure the Snakemake pipeline has run and generated per-gene diagnostic parameters.")
        sys.exit(1)
    
    print("\n=== ANALYSIS COMPLETE ===")
    if results:
        first_result = next(iter(results.values()))
        print(f"Generated separate analyses for {len(results)} POSSE variants")
        print(f"Genes analyzed per variant: {first_result['n_genes']}")
        print(f"Output saved to: {args.output_dir}")
    else:
        print("No analyses completed successfully")
        sys.exit(1)

if __name__ == "__main__":
    main()