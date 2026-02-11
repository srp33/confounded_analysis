#!/usr/bin/env python3
"""
PACE Validation Suite: Simpson's Paradox Stress Test

This suite validates PACE against three distinct failure modes:
1. Over-correction of Biology (Simpson's Paradox)
2. Under-correction of Artifacts (Technical batch effects)
3. Instability (Null hypothesis test)
"""

import numpy as np
import pandas as pd
import sys
import os
from dataclasses import replace
from typing import Dict, Tuple, Optional

# Add the scripts directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pace import PACE_v22, BatchData, PACEHyperparameters

try:
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    print("Warning: sklearn not available, using simple clustering fallback")

def simple_kmeans_fallback(X_pca: np.ndarray, n_clusters: int = 2, random_state: int = 42) -> np.ndarray:
    """Simple K-means fallback when sklearn is not available"""
    np.random.seed(random_state)
    n_samples = X_pca.shape[0]
    
    # Initialize centroids randomly
    centroids = np.random.randn(n_clusters, X_pca.shape[1])
    labels = np.zeros(n_samples, dtype=int)
    
    # Simple K-means iterations
    for _ in range(10):
        # Assign points to nearest centroid
        distances = np.sqrt(((X_pca[:, np.newaxis, :] - centroids[np.newaxis, :, :]) ** 2).sum(axis=2))
        labels = np.argmin(distances, axis=1)
        
        # Update centroids
        for k in range(n_clusters):
            mask = labels == k
            if np.sum(mask) > 0:
                centroids[k] = np.mean(X_pca[mask], axis=0)
    
    return labels

def simple_pca_fallback(X: np.ndarray, n_components: int = 10) -> np.ndarray:
    """Simple PCA fallback when sklearn is not available"""
    # Center the data
    X_centered = X - np.mean(X, axis=0, keepdims=True)
    
    # Compute covariance matrix
    cov_matrix = np.cov(X_centered.T)
    
    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
    
    # Sort by eigenvalues (descending)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvectors = eigenvectors[:, idx]
    
    # Project data
    return X_centered @ eigenvectors[:, :n_components]

def identify_biological_states(X: np.ndarray, n_clusters: int = 2, random_state: int = 42) -> Tuple[np.ndarray, np.ndarray]:
    """
    Identify biological states (cell types) in the data using PCA + clustering
    
    Args:
        X: Gene expression data (genes x samples)
        n_clusters: Number of biological states to identify
        random_state: Random seed for reproducibility
    
    Returns:
        Tuple of (labels, state_indices_list)
    """
    if SKLEARN_AVAILABLE:
        # Use sklearn for PCA and clustering
        pca = PCA(n_components=min(10, X.shape[1]-1), random_state=random_state)
        X_pca = pca.fit_transform(X.T)
        
        kmeans = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10)
        labels = kmeans.fit_predict(X_pca)
    else:
        # Use fallback implementations
        X_pca = simple_pca_fallback(X.T, n_components=min(10, X.shape[1]-1))
        labels = simple_kmeans_fallback(X_pca, n_clusters=n_clusters, random_state=random_state)
    
    # Get indices for each state
    state_indices = []
    for k in range(n_clusters):
        state_k_idxs = np.where(labels == k)[0]
        state_indices.append(state_k_idxs)
    
    return labels, state_indices

def run_validation_suite(full_data: BatchData, pace_model: PACE_v22, 
                        n_genes_check: int = 100, random_state: int = 42) -> Dict[str, float]:
    """
    Run the complete PACE validation suite
    
    Args:
        full_data: Complete dataset to split for testing
        pace_model: Initialized PACE model
        n_genes_check: Number of top variable genes to analyze
        random_state: Random seed for reproducibility
    
    Returns:
        Dictionary of validation scores
    """
    np.random.seed(random_state)
    
    # 0. Pre-requisite: Define Biological States
    X = full_data.data
    labels, state_indices = identify_biological_states(X, n_clusters=2, random_state=random_state)
    
    state_0_idxs, state_1_idxs = state_indices
    
    # Identify top variable genes for metrics
    gene_vars = np.var(X, axis=1)
    top_genes = np.argsort(gene_vars)[-n_genes_check:]
    
    results = {}
    
    # ==========================================
    # Experiment 1: The Biased Split (Simpson's Paradox)
    # ==========================================
    
    # Construct Biased Batches
    n_half = X.shape[1] // 2
    
    # Batch A: 80% State 0, 20% State 1
    n_A_state0 = min(int(n_half * 0.8), len(state_0_idxs))
    n_A_state1 = min(int(n_half * 0.2), len(state_1_idxs))
    
    idx_A = np.concatenate([
        np.random.choice(state_0_idxs, n_A_state0, replace=False),
        np.random.choice(state_1_idxs, n_A_state1, replace=False)
    ])
    
    # Batch B: 20% State 0, 80% State 1
    remaining_state0 = np.setdiff1d(state_0_idxs, idx_A)
    remaining_state1 = np.setdiff1d(state_1_idxs, idx_A)
    
    n_B_state0 = min(int(n_half * 0.2), len(remaining_state0))
    n_B_state1 = min(int(n_half * 0.8), len(remaining_state1))
    
    idx_B = np.concatenate([
        np.random.choice(remaining_state0, n_B_state0, replace=False),
        np.random.choice(remaining_state1, n_B_state1, replace=False)
    ])
    
    # Create batch data
    data_A = BatchData(data=X[:, idx_A], gene_indices=full_data.gene_indices)
    data_B = BatchData(data=X[:, idx_B], gene_indices=full_data.gene_indices)
    
    # Calculate Biological Signal (expected difference due to composition)
    mu_A = np.mean(data_A.data[top_genes], axis=1)
    mu_B = np.mean(data_B.data[top_genes], axis=1)
    bio_diff = mu_B - mu_A
    
    # Run PACE correction
    try:
        corr_B, metadata = pace_model.align(data_A, data_B)
        mu_B_corr = np.mean(corr_B[top_genes], axis=1)
        
        # Metric: Signal Preservation
        # How much of the biological difference survived?
        mask = np.abs(bio_diff) > 0.1  # Filter for genes with real signal
        if np.sum(mask) > 0:
            final_diff = mu_B_corr - mu_A
            preservation = final_diff[mask] / (bio_diff[mask] + 1e-8)  # Avoid division by zero
            score_1 = np.median(preservation)
            
            # Additional metrics
            preservation_mean = np.mean(preservation)
            preservation_std = np.std(preservation)
            n_preserved = np.sum(np.abs(preservation) > 0.8)  # Genes with >80% preservation
            
            results['Exp1_Signal_Preservation_Median'] = score_1
            results['Exp1_Signal_Preservation_Mean'] = preservation_mean
            results['Exp1_Signal_Preservation_Std'] = preservation_std
            results['Exp1_Genes_Well_Preserved'] = n_preserved
            results['Exp1_Total_Genes_Analyzed'] = np.sum(mask)
        else:
            score_1 = 0.0
            results['Exp1_Signal_Preservation_Median'] = 0.0
            results['Exp1_Signal_Preservation_Mean'] = 0.0
            results['Exp1_Signal_Preservation_Std'] = 0.0
            results['Exp1_Genes_Well_Preserved'] = 0
            results['Exp1_Total_Genes_Analyzed'] = 0
            
    except Exception as e:
        results['Exp1_Signal_Preservation_Median'] = -999  # Error code
        results['Exp1_Signal_Preservation_Mean'] = -999
        results['Exp1_Signal_Preservation_Std'] = -999
        results['Exp1_Genes_Well_Preserved'] = 0
        results['Exp1_Total_Genes_Analyzed'] = 0

    # ==========================================
    # Experiment 2: Balanced Split + Forced Artifact
    # ==========================================
    
    # Construct Balanced Batches
    all_indices = np.arange(X.shape[1])
    idx_perm = np.random.permutation(all_indices)
    idx_A_bal = idx_perm[:n_half]
    idx_B_bal = idx_perm[n_half:2*n_half]
    
    data_A_bal = BatchData(data=X[:, idx_A_bal], gene_indices=full_data.gene_indices)
    
    # Add Artificial Batch Effect to B
    artifact_strength = 2.0
    data_B_dirty = BatchData(
        data=X[:, idx_B_bal] + artifact_strength, 
        gene_indices=full_data.gene_indices
    )
    
    # Run PACE correction
    try:
        corr_B_clean, metadata = pace_model.align(data_A_bal, data_B_dirty)
        
        # Metric: Artifact Removal
        # How much of the +2.0 shift was removed?
        original_B = X[:, idx_B_bal]
        correction_applied = np.mean(data_B_dirty.data[top_genes] - corr_B_clean[top_genes], axis=1)
        removal_ratio = correction_applied / artifact_strength
        
        score_2 = np.median(removal_ratio)
        removal_mean = np.mean(removal_ratio)
        removal_std = np.std(removal_ratio)
        n_well_corrected = np.sum(np.abs(removal_ratio - 1.0) < 0.2)  # Within 20% of perfect correction
        
        results['Exp2_Artifact_Removal_Median'] = score_2
        results['Exp2_Artifact_Removal_Mean'] = removal_mean
        results['Exp2_Artifact_Removal_Std'] = removal_std
        results['Exp2_Genes_Well_Corrected'] = n_well_corrected
        results['Exp2_Total_Genes_Analyzed'] = len(top_genes)
        
    except Exception as e:
        results['Exp2_Artifact_Removal_Median'] = -999
        results['Exp2_Artifact_Removal_Mean'] = -999
        results['Exp2_Artifact_Removal_Std'] = -999
        results['Exp2_Genes_Well_Corrected'] = 0
        results['Exp2_Total_Genes_Analyzed'] = len(top_genes)

    # ==========================================
    # Experiment 3: Balanced Null (Stability Test)
    # ==========================================
    
    # Use clean balanced data (no artifact)
    data_B_bal = BatchData(data=X[:, idx_B_bal], gene_indices=full_data.gene_indices)
    
    # Run PACE correction
    try:
        corr_B_null, metadata = pace_model.align(data_A_bal, data_B_bal)
        
        # Metric: Stability
        # How much did we change when no change was needed?
        delta = np.abs(corr_B_null[top_genes] - data_B_bal.data[top_genes])
        instability_mae = np.mean(delta)
        instability_max = np.max(delta)
        instability_std = np.std(delta.flatten())
        
        # Count genes with minimal change
        n_stable_genes = np.sum(np.mean(delta, axis=1) < 0.1)
        
        results['Exp3_Instability_MAE'] = instability_mae
        results['Exp3_Instability_Max'] = instability_max
        results['Exp3_Instability_Std'] = instability_std
        results['Exp3_Stable_Genes'] = n_stable_genes
        results['Exp3_Total_Genes_Analyzed'] = len(top_genes)
        
    except Exception as e:
        results['Exp3_Instability_MAE'] = -999
        results['Exp3_Instability_Max'] = -999
        results['Exp3_Instability_Std'] = -999
        results['Exp3_Stable_Genes'] = 0
        results['Exp3_Total_Genes_Analyzed'] = len(top_genes)
    
    return results

def generate_validation_report(results: Dict[str, float], pace_variant: str = "unknown") -> str:
    """Generate a human-readable validation report"""
    
    report = f"""
# PACE Validation Report: {pace_variant}

## Summary Scores

| Experiment | Metric | Score | Target | Status |
|------------|--------|-------|--------|--------|
| **Simpson's Paradox** | Signal Preservation | {results.get('Exp1_Signal_Preservation_Median', 'N/A'):.3f} | 1.000 | {'✅ PASS' if results.get('Exp1_Signal_Preservation_Median', 0) > 0.8 else '❌ FAIL'} |
| **Artifact Removal** | Correction Ratio | {results.get('Exp2_Artifact_Removal_Median', 'N/A'):.3f} | 1.000 | {'✅ PASS' if results.get('Exp2_Artifact_Removal_Median', 0) > 0.8 else '❌ FAIL'} |
| **Stability Test** | Mean Abs. Change | {results.get('Exp3_Instability_MAE', 'N/A'):.4f} | 0.000 | {'✅ PASS' if results.get('Exp3_Instability_MAE', 999) < 0.1 else '❌ FAIL'} |

## Detailed Results

### Experiment 1: Simpson's Paradox Resistance
- **Signal Preservation (Median):** {results.get('Exp1_Signal_Preservation_Median', 'N/A'):.3f}
- **Signal Preservation (Mean):** {results.get('Exp1_Signal_Preservation_Mean', 'N/A'):.3f} ± {results.get('Exp1_Signal_Preservation_Std', 'N/A'):.3f}
- **Well-Preserved Genes:** {results.get('Exp1_Genes_Well_Preserved', 'N/A')}/{results.get('Exp1_Total_Genes_Analyzed', 'N/A')} ({100*results.get('Exp1_Genes_Well_Preserved', 0)/max(1, results.get('Exp1_Total_Genes_Analyzed', 1)):.1f}%)

### Experiment 2: Technical Artifact Correction
- **Artifact Removal (Median):** {results.get('Exp2_Artifact_Removal_Median', 'N/A'):.3f}
- **Artifact Removal (Mean):** {results.get('Exp2_Artifact_Removal_Mean', 'N/A'):.3f} ± {results.get('Exp2_Artifact_Removal_Std', 'N/A'):.3f}
- **Well-Corrected Genes:** {results.get('Exp2_Genes_Well_Corrected', 'N/A')}/{results.get('Exp2_Total_Genes_Analyzed', 'N/A')} ({100*results.get('Exp2_Genes_Well_Corrected', 0)/max(1, results.get('Exp2_Total_Genes_Analyzed', 1)):.1f}%)

### Experiment 3: Null Stability
- **Mean Absolute Change:** {results.get('Exp3_Instability_MAE', 'N/A'):.4f}
- **Maximum Change:** {results.get('Exp3_Instability_Max', 'N/A'):.4f}
- **Stable Genes:** {results.get('Exp3_Stable_Genes', 'N/A')}/{results.get('Exp3_Total_Genes_Analyzed', 'N/A')} ({100*results.get('Exp3_Stable_Genes', 0)/max(1, results.get('Exp3_Total_Genes_Analyzed', 1)):.1f}%)

## Interpretation

**Signal Preservation (Exp 1):** Measures PACE's ability to resist over-correcting biological differences.
- **1.0** = Perfect preservation of cell-type differences despite compositional bias
- **0.0** = Complete signal loss (like ComBat would do)

**Artifact Removal (Exp 2):** Measures PACE's ability to correct true technical batch effects.
- **1.0** = Perfect removal of artificial +2.0 shift
- **0.0** = No correction applied

**Stability (Exp 3):** Measures PACE's tendency to hallucinate corrections when none are needed.
- **0.0** = Perfect stability (no unnecessary changes)
- **High values** = Unstable, making corrections when none are needed

"""
    
    return report

def main():
    """Main function for command-line usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description="PACE Validation Suite")
    parser.add_argument("--data", required=True, help="Path to input data (RData or CSV)")
    parser.add_argument("--variant", default="default", help="PACE variant to test")
    parser.add_argument("--output", help="Output path for validation report")
    parser.add_argument("--n-genes", type=int, default=100, help="Number of top variable genes to analyze")
    parser.add_argument("--tau", type=float, help="Override tau parameter")
    parser.add_argument("--w-prior", type=float, help="Override w_prior parameter")
    
    args = parser.parse_args()
    
    # Load data (implementation depends on data format)
    print(f"Loading data from {args.data}...")
    # TODO: Implement data loading based on file extension
    
    # Initialize PACE model
    hyperparams = PACEHyperparameters()
    if args.tau is not None:
        hyperparams.tau = args.tau
    if args.w_prior is not None:
        hyperparams.w_prior = args.w_prior
    
    pace_model = PACE_v22(pathway_source='hallmark', hyperparams=hyperparams)
    
    # Run validation
    results = run_validation_suite(full_data, pace_model, n_genes_check=args.n_genes)
    
    # Generate report
    report = generate_validation_report(results, args.variant)
    
    if args.output:
        with open(args.output, 'w') as f:
            f.write(report)
        print(f"Validation report saved to {args.output}")
    else:
        print(report)

if __name__ == "__main__":
    main()