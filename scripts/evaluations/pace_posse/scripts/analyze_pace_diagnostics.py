#!/usr/bin/env python3
"""
Analyze PACE diagnostic outputs to understand performance issues
"""

import os
import re
import pandas as pd
import numpy as np
from pathlib import Path

def extract_pace_metrics_from_logs(log_dir):
    """Extract PACE metrics from log files"""
    
    log_path = Path(log_dir)
    pace_data = []
    
    # Pattern to match PACE log files
    pace_pattern = re.compile(r'pace_(\w+)_(\w+)_(\d+)_(\w+)\.log')
    
    for log_file in log_path.glob('pace_*.log'):
        match = pace_pattern.match(log_file.name)
        if not match:
            continue
            
        variant, classifier, n_datasets, test_study = match.groups()
        
        try:
            with open(log_file, 'r') as f:
                content = f.read()
                
            # Extract metrics using regex
            metrics = {
                'variant': variant,
                'classifier': classifier,
                'n_datasets': int(n_datasets),
                'test_study': test_study,
                'log_file': str(log_file)
            }
            
            # Extract alpha metrics
            alpha_match = re.search(r'Alpha: mean=([\d\.]+), std=([\d\.]+), min=([\d\.]+), max=([\d\.]+)', content)
            if alpha_match:
                metrics.update({
                    'alpha_mean': float(alpha_match.group(1)),
                    'alpha_std': float(alpha_match.group(2)),
                    'alpha_min': float(alpha_match.group(3)),
                    'alpha_max': float(alpha_match.group(4))
                })
            
            # Extract beta metrics
            beta_match = re.search(r'Beta: mean=([\d\.\-]+), std=([\d\.]+)', content)
            if beta_match:
                metrics.update({
                    'beta_mean': float(beta_match.group(1)),
                    'beta_std': float(beta_match.group(2))
                })
            
            # Extract gene and pathway counts
            genes_match = re.search(r'Genes: common=(\d+), unique=(\d+)', content)
            if genes_match:
                metrics.update({
                    'n_common_genes': int(genes_match.group(1)),
                    'n_unique_genes': int(genes_match.group(2))
                })
            
            pathways_match = re.search(r'Pathways used=(\d+), w_prior=([\d\.]+)', content)
            if pathways_match:
                metrics.update({
                    'n_pathways_used': int(pathways_match.group(1)),
                    'w_prior': float(pathways_match.group(2))
                })
            
            # Extract performance metrics
            auc_match = re.search(r'auc\s+([\d\.]+)', content)
            if auc_match:
                metrics['auc'] = float(auc_match.group(1))
                
            acc_match = re.search(r'acc\s+([\d\.]+)', content)
            if acc_match:
                metrics['accuracy'] = float(acc_match.group(1))
            
            # Extract diagnostic assessment
            if 'TOO CONSERVATIVE' in content:
                metrics['diagnostic'] = 'too_conservative'
            elif 'AGGRESSIVE' in content:
                metrics['diagnostic'] = 'aggressive'
            else:
                metrics['diagnostic'] = 'balanced'
                
            pace_data.append(metrics)
            
        except Exception as e:
            print(f"Error processing {log_file}: {e}")
            continue
    
    return pd.DataFrame(pace_data)

def analyze_pace_performance(df):
    """Analyze PACE performance patterns"""
    
    print("=== PACE Performance Analysis ===\n")
    
    # Overall statistics by variant
    print("1. Performance by Variant:")
    if 'auc' in df.columns:
        variant_stats = df.groupby('variant').agg({
            'auc': ['mean', 'std', 'count'],
            'accuracy': ['mean', 'std'] if 'accuracy' in df.columns else 'count'
        }).round(4)
        print(variant_stats)
    
    print("\n2. Alpha Statistics by Variant:")
    if 'alpha_mean' in df.columns:
        alpha_stats = df.groupby('variant').agg({
            'alpha_mean': ['mean', 'std', 'min', 'max'],
            'alpha_std': ['mean', 'std'],
            'prior_strength': 'first'
        }).round(4)
        print(alpha_stats)
    
    print("\n3. Diagnostic Assessment:")
    if 'diagnostic' in df.columns:
        diag_counts = df.groupby(['variant', 'diagnostic']).size().unstack(fill_value=0)
        print(diag_counts)
    
    print("\n4. Pathway Usage:")
    if 'n_pathways_used' in df.columns:
        pathway_stats = df.groupby('variant')['n_pathways_used'].agg(['mean', 'std', 'min', 'max']).round(1)
        print(pathway_stats)
    
    # Identify problematic cases
    print("\n5. Problematic Cases (Alpha mean < 0.15):")
    if 'alpha_mean' in df.columns:
        conservative_cases = df[df['alpha_mean'] < 0.15]
        if not conservative_cases.empty:
            print(f"Found {len(conservative_cases)} cases with very low alpha values:")
            print(conservative_cases[['variant', 'classifier', 'test_study', 'alpha_mean', 'alpha_std', 'auc']].head(10))
        else:
            print("No cases with alpha < 0.15 found")
    
    return df

def generate_recommendations(df):
    """Generate recommendations based on analysis"""
    
    print("\n=== RECOMMENDATIONS ===\n")
    
    if 'alpha_mean' in df.columns:
        # Check if most variants are too conservative
        conservative_fraction = (df['alpha_mean'] < 0.2).mean()
        
        if conservative_fraction > 0.7:
            print("🔴 MAJOR ISSUE: >70% of runs show conservative behavior (alpha < 0.2)")
            print("   Recommendations:")
            print("   - Try even more aggressive settings (prior_strength = 0.0)")
            print("   - Use smaller, more focused pathway sets")
            print("   - Check if pathway-gene overlaps are meaningful")
            print("   - Consider different clamp ranges")
        
        elif conservative_fraction > 0.3:
            print("🟡 MODERATE ISSUE: >30% of runs show conservative behavior")
            print("   Recommendations:")
            print("   - Tune prior_strength parameter")
            print("   - Evaluate pathway relevance to TB biology")
        
        else:
            print("✅ Alpha values appear reasonable across variants")
    
    # Check pathway usage
    if 'n_pathways_used' in df.columns:
        avg_pathways = df['n_pathways_used'].mean()
        if avg_pathways < 10:
            print(f"🔴 LOW PATHWAY USAGE: Average {avg_pathways:.1f} pathways used")
            print("   - Check pathway-gene name matching")
            print("   - Verify gene naming conventions")
        elif avg_pathways > 1000:
            print(f"🟡 HIGH PATHWAY USAGE: Average {avg_pathways:.1f} pathways used")
            print("   - Consider using more focused pathway sets")
            print("   - May be diluting signal with irrelevant pathways")

if __name__ == "__main__":
    # Analyze logs
    log_dir = "outputs/logs/classify_adjusters"
    
    if not os.path.exists(log_dir):
        print(f"Log directory {log_dir} not found")
        exit(1)
    
    print("Extracting PACE metrics from log files...")
    df = extract_pace_metrics_from_logs(log_dir)
    
    if df.empty:
        print("No PACE log files found or no metrics extracted")
        exit(1)
    
    print(f"Found {len(df)} PACE runs to analyze\n")
    
    # Perform analysis
    df = analyze_pace_performance(df)
    
    # Generate recommendations
    generate_recommendations(df)
    
    # Save results
    output_file = "outputs/pace_diagnostic_analysis.csv"
    df.to_csv(output_file, index=False)
    print(f"\nDetailed results saved to: {output_file}")