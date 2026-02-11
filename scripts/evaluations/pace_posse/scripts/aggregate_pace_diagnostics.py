#!/usr/bin/env python3
"""
Aggregate PACE diagnostic information from individual runs.
Generates a summary report of weight sparsity and convergence behavior.
"""

import pandas as pd
import numpy as np
import glob
import os
import sys
import argparse
import json
from pathlib import Path

def load_diagnostic_data(input_dir):
    """Load all diagnostic data from individual runs."""
    
    # Check if input directory exists
    if not os.path.exists(input_dir):
        print(f"ERROR: Diagnostic directory does not exist: {input_dir}", file=sys.stderr)
        return pd.DataFrame()
    
    # Look for CSV files that contain diagnostic information
    pattern = os.path.join(input_dir, "*_diagnostics.csv")
    files = glob.glob(pattern)
    
    print(f"Searching for diagnostic files with pattern: {pattern}")
    print(f"Found {len(files)} diagnostic files")
    
    if not files:
        print(f"No diagnostic files found in {input_dir}")
        print(f"Directory contents:", file=sys.stderr)
        try:
            for item in os.listdir(input_dir):
                print(f"  {item}", file=sys.stderr)
        except Exception as e:
            print(f"  Could not list directory: {e}", file=sys.stderr)
        return pd.DataFrame()
    
    all_data = []
    
    for file in files:
        try:
            print(f"Processing diagnostic file: {os.path.basename(file)}")
            
            # Extract run parameters from filename
            basename = os.path.basename(file).replace('_diagnostics.csv', '')
            parts = basename.split('_')
            
            print(f"  Parsing basename: {basename}")
            print(f"  Parts: {parts}")
            
            if len(parts) >= 4:
                # Handle pace variants like "pace_aggressive_xgboost_4_GSE37250_M"
                if parts[0] == 'pace' and len(parts) >= 5:
                    adjuster = f"{parts[0]}_{parts[1]}"  # pace_aggressive, pace_default, etc.
                    classifier = parts[2]
                    n_datasets = parts[3]
                    test_study = '_'.join(parts[4:])  # Handle multi-part study names
                    print(f"  PACE variant detected: {adjuster}")
                else:
                    adjuster = parts[0]
                    classifier = parts[1] 
                    n_datasets = parts[2]
                    test_study = '_'.join(parts[3:])  # Handle multi-part study names
                    print(f"  Standard adjuster: {adjuster}")
            else:
                # Fallback parsing
                adjuster = "unknown"
                classifier = "unknown"
                n_datasets = "unknown"
                test_study = "unknown"
                print(f"Warning: Could not parse filename {basename}, using defaults")
            
            # Load the diagnostic data
            df = pd.read_csv(file)
            
            if df.empty:
                print(f"Warning: Empty diagnostic file: {file}")
                continue
            
            print(f"  Loaded {len(df)} diagnostic records")
            
            # Add run metadata
            df['adjuster'] = adjuster
            df['classifier'] = classifier
            df['n_datasets'] = n_datasets
            df['test_study'] = test_study
            df['run_id'] = basename
            
            all_data.append(df)
            
        except Exception as e:
            print(f"Error processing {file}: {e}", file=sys.stderr)
            continue
    
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        print(f"Successfully loaded {len(combined_df)} total diagnostic records from {len(all_data)} files")
        return combined_df
    else:
        print("No valid diagnostic data could be loaded from any files")
        return pd.DataFrame()

def generate_diagnostic_summary(df):
    """Generate summary statistics for diagnostic data."""
    
    if df.empty:
        return pd.DataFrame()
    
    # Group by run and calculate final iteration stats
    final_stats = df.groupby('run_id').last().reset_index()
    
    # Calculate summary metrics using CORRECTED interpretation
    summary_stats = []
    
    for adjuster in final_stats['adjuster'].unique():
        adjuster_data = final_stats[final_stats['adjuster'] == adjuster]
        
        stats = {
            'adjuster': adjuster,
            'n_runs': len(adjuster_data),
            # CORRECTED: Use gini_x (reference sparsity) as the main metric
            'avg_final_neff_x_ratio': adjuster_data['neff_x_ratio'].mean(),
            'std_final_neff_x_ratio': adjuster_data['neff_x_ratio'].std(),
            'avg_final_gini_x': adjuster_data['gini_x'].mean(),
            'std_final_gini_x': adjuster_data['gini_x'].std(),
            'avg_iterations': adjuster_data['iter'].mean(),
            # CORRECTED: Low neff_x_ratio is good (specific matching)
            'pct_specific_weights': (adjuster_data['neff_x_ratio'] < 0.5).mean() * 100,
            'pct_sparse_weights': (adjuster_data['gini_x'] > 0.3).mean() * 100
        }
        
        summary_stats.append(stats)
    
    return pd.DataFrame(summary_stats)

def generate_convergence_report(df):
    """Generate convergence behavior analysis."""
    
    if df.empty:
        return ""
    
    report = []
    report.append("# PACE Diagnostic Report")
    report.append("")
    
    # Overall statistics
    total_runs = df['run_id'].nunique()
    report.append(f"**Total Runs Analyzed:** {total_runs}")
    report.append("")
    
    # Final iteration analysis
    final_stats = df.groupby('run_id').last()
    
    # CORRECTED: Weight sparsity analysis using reference sparsity (gini_x)
    specific_pct = (final_stats['neff_x_ratio'] < 0.5).mean() * 100  # Good: using <50% of reference
    sparse_pct = (final_stats['gini_x'] > 0.3).mean() * 100  # Good: reference weights are sparse
    
    report.append("## Reference Matching Analysis")
    report.append(f"- **Specific matching (Neff_X < 50%):** {specific_pct:.1f}% of runs")
    report.append(f"- **Sparse reference weights (Gini_X > 0.3):** {sparse_pct:.1f}% of runs")
    report.append("")
    
    if specific_pct < 70:
        report.append("⚠️  **WARNING:** Low percentage of runs showing specific reference matching!")
        report.append("   This suggests PACE is using too many reference samples rather than matching specific cell types.")
        report.append("")
    else:
        report.append("✅ **GOOD:** Most runs show specific reference matching!")
        report.append("")
    
    # Per-adjuster breakdown
    report.append("## Per-Adjuster Performance")
    report.append("")
    
    for adjuster in sorted(final_stats['adjuster'].unique()):
        adj_data = final_stats[final_stats['adjuster'] == adjuster]
        
        # CORRECTED: Use reference sparsity metrics
        avg_neff_x = adj_data['neff_x_ratio'].mean()
        avg_gini_x = adj_data['gini_x'].mean()
        specific_runs = (adj_data['neff_x_ratio'] < 0.5).sum()
        total_adj_runs = len(adj_data)
        
        report.append(f"### {adjuster}")
        report.append(f"- Average Reference Usage (Neff_X): {avg_neff_x:.1%}")
        report.append(f"- Average Reference Sparsity (Gini_X): {avg_gini_x:.3f}")
        report.append(f"- Specific matching runs: {specific_runs}/{total_adj_runs} ({specific_runs/total_adj_runs*100:.1f}%)")
        
        if avg_neff_x > 0.5:
            report.append("  - ❌ **NEEDS TUNING:** Using too many reference samples")
        elif avg_gini_x > 0.3:
            report.append("  - ✅ **GOOD:** Achieving specific reference matching")
        else:
            report.append("  - ⚠️  **MODERATE:** Some reference specificity detected")
        
        report.append("")
    
    # CORRECTED: Recommendations based on reference sparsity
    report.append("## Recommendations")
    report.append("")
    
    if specific_pct < 70:
        report.append("**Higher Contrast Tuning Needed:**")
        report.append("- Increase tau parameter (e.g., 20.0 → 50.0)")
        report.append("- Reduce w_prior for more aggressive matching")
        report.append("- Verify pathway-based similarity computation")
        report.append("")
    else:
        report.append("**Current settings achieve good reference specificity.**")
        report.append("- Reference matching is working as intended")
        report.append("- Consider fine-tuning tau for optimal performance")
        report.append("")
    
    return "\n".join(report)

def main():
    parser = argparse.ArgumentParser(description='Aggregate PACE diagnostic data')
    parser.add_argument('--input-dir', required=True, help='Directory containing diagnostic CSV files')
    parser.add_argument('--output-summary', required=True, help='Output CSV file for summary statistics')
    parser.add_argument('--output-report', required=True, help='Output markdown file for diagnostic report')
    
    args = parser.parse_args()
    
    # Load diagnostic data
    print(f"Loading diagnostic data from {args.input_dir}")
    df = load_diagnostic_data(args.input_dir)
    
    if df.empty:
        error_msg = f"CRITICAL ERROR: No diagnostic data found in {args.input_dir}!"
        print(error_msg, file=sys.stderr)
        print("Expected diagnostic CSV files with pattern: *_diagnostics.csv", file=sys.stderr)
        print("This indicates PACE runs are not generating diagnostic output.", file=sys.stderr)
        print("Check that:", file=sys.stderr)
        print("  1. PACE wrapper is passing diagnostic_output_path parameter", file=sys.stderr)
        print("  2. PACE.align() method is being called with diagnostic_output_path", file=sys.stderr)
        print("  3. Diagnostic directory has write permissions", file=sys.stderr)
        sys.exit(1)  # Exit with error code to stop Snakemake
    
    print(f"Loaded {len(df)} diagnostic records from {df['run_id'].nunique()} runs")
    
    # Generate summary statistics
    summary_df = generate_diagnostic_summary(df)
    summary_df.to_csv(args.output_summary, index=False)
    print(f"Summary statistics saved to {args.output_summary}")
    
    # Generate diagnostic report
    report = generate_convergence_report(df)
    with open(args.output_report, 'w') as f:
        f.write(report)
    print(f"Diagnostic report saved to {args.output_report}")

if __name__ == "__main__":
    main()