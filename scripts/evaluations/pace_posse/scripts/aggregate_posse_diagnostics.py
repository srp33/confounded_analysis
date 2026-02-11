#!/usr/bin/env python3
"""
Aggregate POSSE v4.0 diagnostic information from individual runs.
Generates a summary report of anchor trust and parameter distributions.
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
                # Handle posse variants like "posse_aggressive_xgboost_4_GSE37250_M"
                if parts[0] == 'posse' and len(parts) >= 5:
                    adjuster = f"{parts[0]}_{parts[1]}"  # posse_aggressive, posse_default, etc.
                    classifier = parts[2]
                    n_datasets = parts[3]
                    test_study = '_'.join(parts[4:])  # Handle multi-part study names
                    print(f"  POSSE variant detected: {adjuster}")
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
    
    # Group by run and calculate summary stats
    summary_stats = []
    
    for adjuster in df['adjuster'].unique():
        adjuster_data = df[df['adjuster'] == adjuster]
        
        # Group by run_id to get per-run statistics
        run_groups = adjuster_data.groupby('run_id')
        
        stats = {
            'adjuster': adjuster,
            'n_runs': len(run_groups),
            'n_genes_total': len(adjuster_data),
            # v4.0 Diagnostic scale factor statistics (with fallback)
            'avg_S_diagnostic': adjuster_data['S_diagnostic'].mean() if 'S_diagnostic' in adjuster_data.columns else np.nan,
            'std_S_diagnostic': adjuster_data['S_diagnostic'].std() if 'S_diagnostic' in adjuster_data.columns else np.nan,
            'min_S_diagnostic': adjuster_data['S_diagnostic'].min() if 'S_diagnostic' in adjuster_data.columns else np.nan,
            'max_S_diagnostic': adjuster_data['S_diagnostic'].max() if 'S_diagnostic' in adjuster_data.columns else np.nan,
            # Trust and alpha statistics from final iteration
            'avg_trust_final': adjuster_data['avg_trust_final'].mean() if 'avg_trust_final' in adjuster_data.columns else np.nan,
            'std_trust_final': adjuster_data['avg_trust_final'].std() if 'avg_trust_final' in adjuster_data.columns else np.nan,
            'avg_alpha_final': adjuster_data['avg_alpha_final'].mean() if 'avg_alpha_final' in adjuster_data.columns else (adjuster_data['alpha_final'].mean() if 'alpha_final' in adjuster_data.columns else np.nan),
            'std_alpha_final': adjuster_data['avg_alpha_final'].std() if 'avg_alpha_final' in adjuster_data.columns else (adjuster_data['alpha_final'].std() if 'alpha_final' in adjuster_data.columns else np.nan),
            # Final parameter statistics - use alpha_final/beta_final as fallback
            'alpha_final_mean': adjuster_data['alpha_final_mean'].mean() if 'alpha_final_mean' in adjuster_data.columns else (adjuster_data['alpha_final'].mean() if 'alpha_final' in adjuster_data.columns else np.nan),
            'beta_final_mean': adjuster_data['beta_final_mean'].mean() if 'beta_final_mean' in adjuster_data.columns else (adjuster_data['beta_final'].mean() if 'beta_final' in adjuster_data.columns else np.nan),
            # Configuration parameters (simplified v4.0)
            'avg_tau': adjuster_data['tau'].mean() if 'tau' in adjuster_data.columns else 20.0,
            'avg_top_k_percent': adjuster_data['top_k_percent'].mean() if 'top_k_percent' in adjuster_data.columns else 0.15,
            # Correction assessment - use alpha_final as fallback
            'correction_magnitude': abs(adjuster_data['alpha_final_mean'].mean() - 1.0) if 'alpha_final_mean' in adjuster_data.columns else (abs(adjuster_data['alpha_final'].mean() - 1.0) if 'alpha_final' in adjuster_data.columns else np.nan),
            'shift_magnitude': abs(adjuster_data['beta_final_mean'].mean()) if 'beta_final_mean' in adjuster_data.columns else (abs(adjuster_data['beta_final'].mean()) if 'beta_final' in adjuster_data.columns else np.nan),
            # Gene coverage (v4.0 specific)
            'avg_genes_covered': adjuster_data['genes_covered'].mean() if 'genes_covered' in adjuster_data.columns else np.nan,
            'avg_genes_total': adjuster_data['genes_total'].mean() if 'genes_total' in adjuster_data.columns else np.nan
        }
        
        summary_stats.append(stats)
    
    return pd.DataFrame(summary_stats)

def generate_convergence_report(df):
    """Generate POSSE v4.1 diagnostic analysis report with timidity analysis."""
    
    if df.empty:
        return ""
    
    report = []
    report.append("# POSSE v4.1 Diagnostic Report")
    report.append("## Robust HK + Pure Local Trust Analysis")
    report.append("")
    
    # Overall statistics
    total_runs = df['run_id'].nunique()
    total_genes = len(df)
    report.append(f"**Total Runs Analyzed:** {total_runs}")
    report.append(f"**Total Gene Records:** {total_genes}")
    report.append("")
    
    # Check available columns and adapt
    has_s_diagnostic = 'S_diagnostic' in df.columns
    has_alpha_final_mean = 'alpha_final_mean' in df.columns
    has_alpha_final = 'alpha_final' in df.columns
    
    # Use alpha_final if alpha_final_mean not available
    alpha_col = 'alpha_final_mean' if has_alpha_final_mean else 'alpha_final'
    beta_col = 'beta_final_mean' if 'beta_final_mean' in df.columns else 'beta_final'
    
    # Diagnostic scale analysis (v4.1 with HK scaler)
    report.append("## Scale Analysis (Safety Valve + HK Scaler)")
    
    if has_s_diagnostic:
        avg_s_diag = df['S_diagnostic'].mean()
        std_s_diag = df['S_diagnostic'].std()
    else:
        # Estimate scale from alpha_final (inverse relationship)
        avg_s_diag = 1.0 / df[alpha_col].mean() if df[alpha_col].mean() != 0 else 1.0
        std_s_diag = df[alpha_col].std() if has_alpha_final else 0.0
        report.append("*Note: S_diagnostic not available, estimating from alpha_final*")
        report.append("")
    
    report.append(f"- **Average S_diagnostic (HK Scale Factor):** {avg_s_diag:.3f} ± {std_s_diag:.3f}")
    
    if avg_s_diag < 0.01 or avg_s_diag > 100.0:
        report.append("  - ❌ **CRITICAL:** Extreme scale differences detected!")
        report.append("    Safety valve should have triggered CPM normalization.")
        report.append("    Check if library size differences exceed 10x threshold.")
    elif avg_s_diag < 0.1 or avg_s_diag > 10.0:
        report.append("  - ⚠️  **WARNING:** Large technical gains detected!")
        report.append("    HK scaler is handling moderate technical artifacts (2x-10x).")
        report.append("    This is normal for cross-platform or cross-protocol data.")
    elif 0.8 <= avg_s_diag <= 1.2:
        report.append("  - ✅ **EXCELLENT:** HK genes show minimal technical differences.")
        report.append("    Batches are well-matched after safety valve normalization.")
    else:
        report.append("  - ✅ **GOOD:** HK scaler handling moderate technical differences.")
    
    report.append("")
    
    # Correction Effectiveness Analysis
    report.append("## Correction Effectiveness Analysis")
    
    avg_alpha_final = df[alpha_col].mean()
    avg_beta_final = df[beta_col].mean()
    
    report.append(f"- **Final Alpha (Scale Correction):** {avg_alpha_final:.3f}")
    report.append(f"- **Final Beta (Shift Correction):** {avg_beta_final:.3f}")
    
    # Timidity Analysis - Key diagnostic for v4.1
    report.append("")
    report.append("## Timidity Analysis")
    report.append("*Comparing POSSE aggressiveness vs expected ComBat-level corrections*")
    report.append("")
    
    # Calculate correction magnitude
    if 'correction_magnitude' in df.columns:
        avg_correction = df['correction_magnitude'].mean()
        report.append(f"- **POSSE Correction Magnitude:** {avg_correction:.4f}")
        
        # Estimate what ComBat would do (rough heuristic)
        # ComBat typically applies corrections proportional to detected scale differences
        expected_combat_magnitude = abs(1.0 - avg_s_diag) * 0.8  # ComBat ~80% correction
        timidity_ratio = avg_correction / (expected_combat_magnitude + 1e-9)
        
        report.append(f"- **Expected ComBat Magnitude:** {expected_combat_magnitude:.4f}")
        report.append(f"- **Timidity Ratio (POSSE/ComBat):** {timidity_ratio:.3f}")
        
        if timidity_ratio < 0.3:
            report.append("  - ❌ **SEVERELY TIMID:** POSSE applying <30% of expected corrections")
            report.append("    - Trust gating is blocking necessary corrections")
            report.append("    - ComBat likely outperforming due to aggressive approach")
            report.append("    - **Action:** Increase tau parameter for higher confidence")
        elif timidity_ratio < 0.6:
            report.append("  - ⚠️  **MODERATELY TIMID:** POSSE applying <60% of expected corrections")
            report.append("    - Algorithm being conservative with biological preservation")
            report.append("    - May lose to ComBat on technical noise removal")
        elif timidity_ratio > 1.5:
            report.append("  - ⚠️  **OVER-AGGRESSIVE:** POSSE applying >150% of expected corrections")
            report.append("    - Risk of removing biological signal")
            report.append("    - Consider reducing tau parameter")
        else:
            report.append("  - ✅ **BALANCED:** POSSE correction magnitude appropriate")
            report.append("    - Good balance between artifact removal and biology preservation")
    
    # Check if algorithm is actually correcting the detected scale difference
    expected_alpha = 1.0 / avg_s_diag  # What alpha should be to correct S_diagnostic
    alpha_error = abs(avg_alpha_final - expected_alpha)
    
    report.append(f"- **Expected Alpha (1/S_diagnostic):** {expected_alpha:.3f}")
    report.append(f"- **Alpha Error:** {alpha_error:.3f}")
    
    if alpha_error < 0.1:
        report.append("  - ✅ **SUCCESS:** HK scaler + local trust correctly removing scale artifact")
    elif alpha_error < 0.3:
        report.append("  - ⚠️  **PARTIAL:** Algorithm partially correcting the scale artifact")
        report.append("    - Trust gating may be dampening some corrections")
    else:
        report.append("  - ❌ **FAILURE:** Algorithm not correcting the detected scale artifact")
        report.append("    - HK scaler detected the problem but local trust blocked correction")
        report.append("    - **Action:** Check pathway gene coverage and correlation quality")
    
    report.append("")
    
    # Per-adjuster breakdown
    report.append("## Per-Adjuster Performance")
    report.append("")
    
    for adjuster in sorted(df['adjuster'].unique()):
        adj_data = df[df['adjuster'] == adjuster]
        
        if has_s_diagnostic:
            adj_s_diag = adj_data['S_diagnostic'].mean()
        else:
            adj_s_diag = 1.0 / adj_data[alpha_col].mean() if adj_data[alpha_col].mean() != 0 else 1.0
        adj_alpha = adj_data[alpha_col].mean()
        adj_beta = adj_data[beta_col].mean()
        
        report.append(f"### {adjuster}")
        report.append(f"- HK Scale Factor (S_diagnostic): {adj_s_diag:.3f}")
        report.append(f"- Final Alpha: {adj_alpha:.3f}")
        report.append(f"- Final Beta: {adj_beta:.3f}")
        
        if 'correction_magnitude' in adj_data.columns:
            adj_correction = adj_data['correction_magnitude'].mean()
            report.append(f"- Correction Magnitude: {adj_correction:.4f}")
            
            # Variant-specific timidity analysis
            expected_combat = abs(1.0 - adj_s_diag) * 0.8
            variant_timidity = adj_correction / (expected_combat + 1e-9)
            report.append(f"- Timidity Ratio: {variant_timidity:.3f}")
        
        # Performance assessment for this adjuster
        expected_alpha = 1.0 / adj_s_diag
        alpha_error = abs(adj_alpha - expected_alpha)
        
        if alpha_error < 0.2:
            report.append("  - ✅ **EXCELLENT:** Accurate HK-based correction")
        elif alpha_error > 0.5:
            report.append("  - ❌ **NEEDS ATTENTION:** Poor correction accuracy")
            if 'correction_magnitude' in adj_data.columns and variant_timidity < 0.4:
                report.append("    - Likely cause: Excessive timidity (trust gating too strict)")
        else:
            report.append("  - ✅ **GOOD:** Reasonable correction performance")
        
        report.append("")
    
    # Configuration analysis
    report.append("## Configuration Analysis")
    report.append("")
    
    # Get configuration values for v4.1
    config_cols = ['tau', 'top_k_percent', 'pathway_source']
    available_cols = [col for col in config_cols if col in df.columns]
    
    if available_cols:
        unique_configs = df[['adjuster'] + available_cols].drop_duplicates()
        
        for _, config in unique_configs.iterrows():
            report.append(f"### {config['adjuster']} Configuration")
            for col in available_cols:
                if col in config:
                    report.append(f"- {col}: {config[col]}")
            report.append("")
    
    # Enhanced Recommendations for v4.1
    report.append("## Recommendations")
    report.append("")
    
    if avg_s_diag < 0.01 or avg_s_diag > 100.0:
        report.append("**CRITICAL Library Size Mismatch:**")
        report.append("- Safety valve should trigger at 10x differences")
        report.append("- Verify CPM normalization is working correctly")
        report.append("- Check for raw counts vs normalized data mixing")
        report.append("")
    elif 'correction_magnitude' in df.columns:
        avg_correction = df['correction_magnitude'].mean()
        expected_combat = abs(1.0 - avg_s_diag) * 0.8
        timidity_ratio = avg_correction / (expected_combat + 1e-9)
        
        if timidity_ratio < 0.3:
            report.append("**Severe Timidity - Algorithm Too Conservative:**")
            report.append("- POSSE is under-correcting compared to ComBat")
            report.append("- **Primary Action:** Increase tau parameter (try 30.0-50.0)")
            report.append("- **Secondary:** Reduce top_k_percent for more selective peers")
            report.append("- **Verify:** Pathway gene coverage matches your data")
            report.append("")
        elif timidity_ratio < 0.6:
            report.append("**Moderate Timidity - Room for Improvement:**")
            report.append("- POSSE being cautious, may lose to ComBat on technical noise")
            report.append("- **Try:** Increase tau parameter (try 25.0-35.0)")
            report.append("- **Monitor:** Signal preservation vs artifact removal balance")
            report.append("")
        elif timidity_ratio > 1.5:
            report.append("**Over-Aggressive Correction:**")
            report.append("- Risk of removing biological signal")
            report.append("- **Action:** Reduce tau parameter (try 10.0-15.0)")
            report.append("- **Monitor:** Ensure biological differences are preserved")
            report.append("")
        else:
            report.append("**Optimal Performance Achieved:**")
            report.append("- POSSE v4.1 balancing artifact removal with biology preservation")
            report.append("- HK scaler handling technical gains appropriately")
            report.append("- Local trust providing surgical precision")
            report.append("- Current parameters appear well-tuned")
            report.append("")
    
    # v4.1 specific guidance
    report.append("## POSSE v4.1 Specific Notes")
    report.append("")
    report.append("**Safety Valve (CPM Normalization):**")
    report.append("- Triggers automatically for >10x library size differences")
    report.append("- Handles raw counts vs normalized data mismatches")
    report.append("")
    report.append("**Robust HK Scaler:**")
    report.append("- Uses bottom 20% CoV genes as housekeeping anchors")
    report.append("- Distinguishes technical gains from biological bias")
    report.append("- Pre-scales data before local peer finding")
    report.append("")
    report.append("**Pure Local Trust:**")
    report.append("- Max-pooling of correlation and stability scores")
    report.append("- No magic number parameters (anchor_trust_bonus removed)")
    report.append("- Surgical correction with biological preservation")
    report.append("")
    
    return "\n".join(report)

def main():
    parser = argparse.ArgumentParser(description='Aggregate POSSE v4.0 diagnostic data')
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
        print("This indicates POSSE runs are not generating diagnostic output.", file=sys.stderr)
        print("Check that:", file=sys.stderr)
        print("  1. POSSE wrapper is passing save_diagnostics parameter", file=sys.stderr)
        print("  2. POSSE.align() method is generating diagnostic output", file=sys.stderr)
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