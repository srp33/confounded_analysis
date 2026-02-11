#!/usr/bin/env python3
"""
Pathway Coherence Analysis for PACE validation.
Validates the core intuition that genes within pathways move together.
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import sys
import os
import argparse
from pathlib import Path

# Add the scripts directory to path to import pace
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pace import load_pathways, BatchData

def analyze_pathway_coherence(data_obj, pathway_dict, min_genes=5):
    """
    Validates the 'Core Intuition' by measuring how well each gene 
    correlates with the consensus (Eigengene) of its parent pathway.
    
    Args:
        data_obj: BatchData object (contains .data and .gene_indices)
        pathway_dict: Dictionary of {pathway_name: [gene_names]}
        min_genes: Minimum genes required for pathway analysis
    
    Returns:
        df_coherence: DataFrame with [Gene, Pathway, Correlation, Variance_Explained]
    """
    print("Starting Pathway Coherence Analysis...")
    
    # Map genes to matrix indices
    gene_map = {g: i for i, g in enumerate(data_obj.gene_indices)}
    results = []
    
    processed_pathways = 0
    skipped_pathways = 0
    
    for pathway, genes in pathway_dict.items():
        # 1. Extract Pathway Matrix
        idxs = [gene_map[g] for g in genes if g in gene_map]
        
        if len(idxs) < min_genes:
            skipped_pathways += 1
            continue
            
        # Matrix X_p: (Genes x Samples)
        X_p = data_obj.data[idxs, :]
        
        # 2. Calculate Eigengene (PC1 of the pathway)
        # We perform PCA on the transpose (Samples x Genes) to get sample embeddings
        # Standardize first (Mean centering is crucial for PCA)
        X_p_centered = X_p - np.mean(X_p, axis=1, keepdims=True)
        
        try:
            # We only need the first component
            pca = PCA(n_components=1)
            # Fit on (Samples x Genes)
            pca.fit(X_p_centered.T)
            
            # The Eigengene is the projection of samples onto PC1
            # Shape: (Samples,)
            eigengene = pca.transform(X_p_centered.T).flatten()
            
            # Flip sign if needed to align with mean expression 
            # (PCA sign is arbitrary, we want positive correlation with activity)
            if np.mean(eigengene * np.mean(X_p_centered, axis=0)) < 0:
                eigengene = -eigengene
                
            # 3. Correlate individual genes with the Eigengene
            # This answers: "Does this gene follow the pathway leader?"
            pathway_variance_explained = pca.explained_variance_ratio_[0]
            
            for i, gene_idx in enumerate(idxs):
                gene_vec = X_p_centered[i, :]
                gene_name = data_obj.gene_indices[gene_idx]
                
                # Pearson Correlation
                if np.std(gene_vec) > 1e-8 and np.std(eigengene) > 1e-8:
                    corr = np.corrcoef(gene_vec, eigengene)[0, 1]
                else:
                    corr = 0.0
                
                results.append({
                    'Gene': gene_name,
                    'Pathway': pathway,
                    'Gene_Pathway_Corr': corr, # How much this gene follows the pathway
                    'Pathway_PC1_Var': pathway_variance_explained, # How coherent the pathway is globally
                    'Pathway_Size': len(idxs)
                })
            
            processed_pathways += 1
            if processed_pathways % 100 == 0:
                print(f"  Processed {processed_pathways} pathways...")
                
        except Exception as e:
            print(f"Skipping {pathway}: {e}")
            skipped_pathways += 1
            continue

    df = pd.DataFrame(results)
    
    # 4. Summary Statistics
    print(f"\n=== PATHWAY COHERENCE ANALYSIS RESULTS ===")
    print(f"Processed pathways: {processed_pathways}")
    print(f"Skipped pathways: {skipped_pathways}")
    print(f"Analyzed {len(df)} Gene-Pathway pairs.")
    
    if len(df) > 0:
        avg_corr = df['Gene_Pathway_Corr'].mean()
        median_corr = df['Gene_Pathway_Corr'].median()
        high_coherence_pct = (df['Gene_Pathway_Corr'] > 0.4).mean() * 100
        low_coherence_pct = (df['Gene_Pathway_Corr'] < 0.1).mean() * 100
        
        print(f"Average Gene-Pathway Coherence: {avg_corr:.3f}")
        print(f"Median Gene-Pathway Coherence: {median_corr:.3f}")
        print(f"High coherence (>0.4): {high_coherence_pct:.1f}%")
        print(f"Low coherence (<0.1): {low_coherence_pct:.1f}%")
        
        # Pathway-level statistics
        pathway_stats = df.groupby('Pathway').agg({
            'Gene_Pathway_Corr': ['mean', 'std', 'count'],
            'Pathway_PC1_Var': 'first'
        }).round(3)
        
        pathway_stats.columns = ['Mean_Corr', 'Std_Corr', 'N_Genes', 'PC1_Var']
        pathway_stats = pathway_stats.sort_values('Mean_Corr', ascending=False)
        
        print(f"\nTop 10 Most Coherent Pathways:")
        print(pathway_stats.head(10))
        
        print(f"\nBottom 10 Least Coherent Pathways:")
        print(pathway_stats.tail(10))
        
        # Interpretation
        print(f"\n=== INTERPRETATION ===")
        if avg_corr > 0.4:
            print("✅ EXCELLENT: High average coherence validates the core PACE intuition!")
            print("   Genes within pathways move together. PACE aggressive mode is safe with proper regularization.")
        elif avg_corr > 0.2:
            print("✅ GOOD: Moderate coherence supports the PACE approach.")
            print("   Most pathways show meaningful structure. PACE should work well.")
        elif avg_corr > 0.1:
            print("⚠️  MODERATE: Low-moderate coherence suggests mixed pathway quality.")
            print("   Some pathways work well, others may be noisy. Consider pathway filtering.")
        else:
            print("❌ POOR: Very low coherence suggests pathway-gene mismatch.")
            print("   Current gene sets may not match your tissue type. Consider different pathway database.")
    
    return df

def main():
    parser = argparse.ArgumentParser(description='Analyze pathway coherence for PACE validation')
    parser.add_argument('--data-file', required=True, help='Path to RData file containing reference data')
    parser.add_argument('--pathway-file', required=True, help='Path to pathway GMT file')
    parser.add_argument('--output-csv', required=True, help='Output CSV file for detailed results')
    parser.add_argument('--output-report', required=True, help='Output text report file')
    parser.add_argument('--min-genes', type=int, default=5, help='Minimum genes per pathway')
    
    args = parser.parse_args()
    
    print("Loading pathway coherence analysis...")
    
    # Load pathways
    print(f"Loading pathways from {args.pathway_file}")
    pathway_dict = load_pathways(args.pathway_file, organism='Human')
    print(f"Loaded {len(pathway_dict)} pathways")
    
    # Load reference data (this would need to be adapted based on your data format)
    print(f"Loading reference data from {args.data_file}")
    # This is a placeholder - you'd need to implement data loading based on your format
    # For now, assume we have a way to load the reference batch
    
    print("ERROR: Data loading not implemented yet.")
    print("This script needs to be integrated with your specific data loading pipeline.")
    print("The analyze_pathway_coherence function is ready to use once you have:")
    print("  1. A BatchData object with reference samples")
    print("  2. The pathway dictionary")
    
    return 1

if __name__ == "__main__":
    main()