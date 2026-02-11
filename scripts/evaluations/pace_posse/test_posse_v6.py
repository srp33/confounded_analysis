#!/usr/bin/env python3
"""Quick test of POSSE v6.0 fixes"""

import sys
sys.path.append('scripts')

import numpy as np
import subprocess
import os

from posse_v6 import POSSEv6, POSSEv6Hyperparameters, BatchData

def load_one_study():
    r_script = """
    load('data/TB_real_data.RData')
    write.table(dat_lst[['Africa']], 'temp_africa.txt', sep='\\t', row.names=FALSE, col.names=FALSE)
    write.table(rownames(dat_lst[['Africa']]), 'temp_genes.txt', row.names=FALSE, col.names=FALSE, quote=FALSE)
    """
    subprocess.run(['R', '--slave', '--vanilla'], input=r_script, text=True, capture_output=True)
    data = np.loadtxt('temp_africa.txt', delimiter='\t')
    genes = np.loadtxt('temp_genes.txt', dtype=str)
    os.remove('temp_africa.txt')
    os.remove('temp_genes.txt')
    return data, genes

def get_pathways(genes):
    pathways = {
        'housekeeping': ['ACTB', 'GAPDH', 'B2M', 'HPRT1', 'TBP', 'YWHAZ'],
        'interferon': ['IFNG', 'IRF1', 'IRF7', 'STAT1', 'GBP1', 'GBP2']
    }
    filtered = {}
    for name, gene_list in pathways.items():
        present = [g for g in gene_list if g in genes]
        if len(present) >= 3:
            filtered[name] = present
    return filtered

def test_silence():
    """Identical data should give α≈1, β≈0"""
    print("TEST 1: SILENCE (identical data)")
    data, genes = load_one_study()
    pathways = get_pathways(genes)
    
    ref = BatchData(data=data, gene_indices=genes)
    target = BatchData(data=data.copy(), gene_indices=genes)
    
    posse = POSSEv6(pathway_dict=pathways)
    corrected, meta = posse.align(ref, target)
    
    rmse = np.sqrt(np.mean((corrected.data - data)**2))
    print(f"  RMSE: {rmse:.6f} (should be ~0)")
    print(f"  Alpha: {meta['alpha_mean']:.4f} (should be ~1)")
    print(f"  Beta: {meta['beta_mean']:.4f} (should be ~0)")
    return rmse < 0.1

def test_shift():
    """Known shift should be recovered"""
    print("\nTEST 2: SHIFT (Y = X + 2.0)")
    data, genes = load_one_study()
    pathways = get_pathways(genes)
    
    shift = 2.0
    shifted = data + shift
    
    ref = BatchData(data=data, gene_indices=genes)
    target = BatchData(data=shifted, gene_indices=genes)
    
    posse = POSSEv6(pathway_dict=pathways)
    corrected, meta = posse.align(ref, target)
    
    rmse = np.sqrt(np.mean((corrected.data - data)**2))
    print(f"  Recovery RMSE: {rmse:.4f} (should be ~0)")
    print(f"  Alpha: {meta['alpha_mean']:.4f} (should be ~1)")
    print(f"  Beta: {meta['beta_mean']:.4f}")
    return rmse < 0.5

def test_scale():
    """Known scale should be recovered"""
    print("\nTEST 3: SCALE (Y = 2.0 * X)")
    data, genes = load_one_study()
    pathways = get_pathways(genes)
    
    scale = 2.0
    scaled = data * scale
    
    ref = BatchData(data=data, gene_indices=genes)
    target = BatchData(data=scaled, gene_indices=genes)
    
    posse = POSSEv6(pathway_dict=pathways)
    corrected, meta = posse.align(ref, target)
    
    rmse = np.sqrt(np.mean((corrected.data - data)**2))
    print(f"  Recovery RMSE: {rmse:.4f} (should be ~0)")
    print(f"  Alpha: {meta['alpha_mean']:.4f} (should be ~0.5)")
    print(f"  Beta: {meta['beta_mean']:.4f}")
    return rmse < 0.5

def test_gene_specific():
    """Gene-specific batch effects"""
    print("\nTEST 4: GENE-SPECIFIC (realistic)")
    data, genes = load_one_study()
    pathways = get_pathways(genes)
    
    n_genes = data.shape[0]
    np.random.seed(42)
    gene_alpha = np.random.normal(1.2, 0.2, size=(n_genes, 1))
    gene_beta = np.random.normal(0.5, 0.5, size=(n_genes, 1))
    
    distorted = data * gene_alpha + gene_beta
    
    ref = BatchData(data=data, gene_indices=genes)
    target = BatchData(data=distorted, gene_indices=genes)
    
    posse = POSSEv6(pathway_dict=pathways)
    corrected, meta = posse.align(ref, target)
    
    rmse = np.sqrt(np.mean((corrected.data - data)**2))
    print(f"  Recovery RMSE: {rmse:.4f}")
    print(f"  Alpha: {meta['alpha_mean']:.4f}")
    print(f"  Beta: {meta['beta_mean']:.4f}")
    return rmse < 1.0

if __name__ == "__main__":
    print("POSSE v6.0 VALIDATION")
    print("=" * 50)
    
    results = {
        'silence': test_silence(),
        'shift': test_shift(),
        'scale': test_scale(),
        'gene_specific': test_gene_specific()
    }
    
    print("\n" + "=" * 50)
    print("SUMMARY")
    for test, passed in results.items():
        print(f"  {test}: {'PASS' if passed else 'FAIL'}")
