#!/usr/bin/env python3
"""
Diagnostic script to visualize the "Timidity Gap" between ComBat and POSSE.
Compares correction aggressiveness and agreement between the two methods.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def compare_adjusters(raw_data, combat_data, posse_data, gene_names=None):
    """
    Visualizes the "Timidity Gap" between ComBat (Aggressive) and POSSE (Conservative).
    
    Args:
        raw_data: Original uncorrected data (genes x samples)
        combat_data: ComBat corrected data (genes x samples)
        posse_data: POSSE corrected data (genes x samples)
        gene_names: Optional gene names for labeling
    """
    # 1. Calculate Absolute Corrections
    diff_combat = np.abs(combat_data - raw_data)
    diff_posse = np.abs(posse_data - raw_data)
    
    mean_diff_c = np.mean(diff_combat, axis=1)
    mean_diff_p = np.mean(diff_posse, axis=1)
    
    # 2. Plot 1: Histogram of Correction Strength
    plt.figure(figsize=(15, 5))
    
    plt.subplot(1, 3, 1)
    plt.hist(mean_diff_c, bins=50, alpha=0.7, label='ComBat', color='red', density=True)
    plt.hist(mean_diff_p, bins=50, alpha=0.7, label='POSSE', color='blue', density=True)
    plt.title("Correction Aggressiveness\n(Mean Absolute Adjustment per Gene)")
    plt.xlabel("Magnitude of Change")
    plt.legend()
    
    # 3. Plot 2: Scatter of Corrections (Do they agree on direction?)
    # We take the mean shift per gene as a proxy for Beta
    beta_c = np.mean(combat_data - raw_data, axis=1)
    beta_p = np.mean(posse_data - raw_data, axis=1)
    
    plt.subplot(1, 3, 2)
    plt.scatter(beta_c, beta_p, alpha=0.5, s=10)
    
    # Add y=x line (Perfect Agreement)
    lims = [
        np.min([plt.xlim(), plt.ylim()]),  
        np.max([plt.xlim(), plt.ylim()]),  
    ]
    plt.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
    plt.title(f"Correction Agreement\ncorr = {np.corrcoef(beta_c, beta_p)[0,1]:.3f}")
    plt.xlabel("ComBat Shift (Beta)")
    plt.ylabel("POSSE Shift (Beta)")
    
    # 4. Plot 3: Ratio of Corrections
    # If Ratio < 1.0, POSSE is timid.
    ratio = (np.abs(beta_p) + 1e-9) / (np.abs(beta_c) + 1e-9)
    
    plt.subplot(1, 3, 3)
    plt.hist(ratio, bins=50, range=(0, 2), color='purple', alpha=0.7)
    plt.axvline(1.0, color='k', linestyle='--')
    plt.title(f"POSSE / ComBat Ratio\nMedian: {np.median(ratio):.2f}")
    plt.xlabel("Ratio (<1.0 means POSSE is timid)")
    
    plt.tight_layout()
    plt.savefig("posse_vs_combat_debug.png")
    print("Diagnostic plot saved to 'posse_vs_combat_debug.png'")
    
    # Print summary statistics
    print("\n=== TIMIDITY ANALYSIS ===")
    print(f"ComBat mean correction magnitude: {np.mean(mean_diff_c):.4f}")
    print(f"POSSE mean correction magnitude: {np.mean(mean_diff_p):.4f}")
    print(f"Correction correlation: {np.corrcoef(beta_c, beta_p)[0,1]:.3f}")
    print(f"Median timidity ratio: {np.median(ratio):.3f}")
    print(f"Fraction of genes where POSSE < ComBat: {np.mean(ratio < 1.0):.3f}")
    
    return {
        'combat_magnitude': mean_diff_c,
        'posse_magnitude': mean_diff_p,
        'combat_shift': beta_c,
        'posse_shift': beta_p,
        'timidity_ratio': ratio,
        'correlation': np.corrcoef(beta_c, beta_p)[0,1],
        'median_ratio': np.median(ratio)
    }

if __name__ == "__main__":
    # Example usage with synthetic data
    print("Testing compare_correctors with synthetic data...")
    
    # Generate test data
    np.random.seed(42)
    n_genes, n_samples = 1000, 100
    
    # Raw data with batch effect
    raw_data = np.random.randn(n_genes, n_samples)
    batch_effect = np.random.randn(n_genes, 1) * 0.5
    raw_data[:, 50:] += batch_effect  # Add batch effect to second half
    
    # Simulate ComBat (aggressive correction)
    combat_data = raw_data.copy()
    combat_data[:, 50:] -= batch_effect * 0.8  # Strong correction
    
    # Simulate POSSE (conservative correction)
    posse_data = raw_data.copy()
    posse_data[:, 50:] -= batch_effect * 0.3  # Weak correction
    
    # Run comparison
    results = compare_adjusters(raw_data, combat_data, posse_data)
    print(f"Test completed. Median timidity ratio: {results['median_ratio']:.3f}")