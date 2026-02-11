#!/usr/bin/env python3
"""Debug POSSEv6 on cross-study data"""

import sys
sys.path.append('scripts')

import numpy as np
import subprocess
import os

from posse_v6 import POSSEv6, POSSEv6Hyperparameters, BatchData as BatchDataV6

def load_two_studies():
    """Load two studies from R"""
    r_script = """
    load('data/TB_real_data.RData')
    write.table(dat_lst[['Africa']], 'temp_train_data.txt', sep='\\t', row.names=FALSE, col.names=FALSE)
    write.table(dat_lst[['GSE37250_M']], 'temp_test_data.txt', sep='\\t', row.names=FALSE, col.names=FALSE)
    write.table(rownames(dat_lst[['Africa']]), 'temp_genes.txt', row.names=FALSE, col.names=FALSE, quote=FALSE)
    """
    subprocess.run(['R', '--slave', '--vanilla'], input=r_script, text=True, capture_output=True)
    
    train_data = np.loadtxt('temp_train_data.txt', delimiter='\t')
    test_data = np.loadtxt('temp_test_data.txt', delimiter='\t')
    genes = np.loadtxt('temp_genes.txt', dtype=str)
    
    for f in ['temp_train_data.txt', 'temp_test_data.txt', 'temp_genes.txt']:
        if os.path.exists(f):
            os.remove(f)
    
    return train_data, test_data, genes

def get_pathways(genes):
    pathways = {
        'housekeeping': ['ACTB', 'GAPDH', 'B2M', 'HPRT1', 'TBP', 'YWHAZ'],
        'interferon': ['IFNG', 'IRF1', 'IRF7', 'STAT1', 'GBP1', 'GBP2'],
        'inflammatory': ['TNF', 'IL1B', 'IL6', 'NFKB1']
    }
    filtered = {}
    for name, gene_list in pathways.items():
        present = [g for g in gene_list if g in genes]
        if len(present) >= 3:
            filtered[name] = present
    return filtered

def main():
    print("DEBUG: POSSEv6 on cross-study data")
    print("=" * 60)
    
    train_data, test_data, genes = load_two_studies()
    pathways = get_pathways(genes)
    
    print(f"Train (Africa): shape={train_data.shape}")
    print(f"  Mean: {np.mean(train_data):.3f}, Std: {np.std(train_data):.3f}")
    print(f"  Range: [{np.min(train_data):.3f}, {np.max(train_data):.3f}]")
    
    print(f"\nTest (GSE37250_M): shape={test_data.shape}")
    print(f"  Mean: {np.mean(test_data):.3f}, Std: {np.std(test_data):.3f}")
    print(f"  Range: [{np.min(test_data):.3f}, {np.max(test_data):.3f}]")
    
    # The data is VERY different - this is the batch effect!
    print(f"\nBatch effect magnitude:")
    print(f"  Mean shift: {np.mean(test_data) - np.mean(train_data):.3f}")
    print(f"  Std ratio: {np.std(test_data) / np.std(train_data):.3f}")
    
    # Run POSSEv6
    ref = BatchDataV6(data=train_data, gene_indices=genes)
    target = BatchDataV6(data=test_data, gene_indices=genes)
    
    posse = POSSEv6(pathway_dict=pathways, hyperparams=POSSEv6Hyperparameters(tau=25.0))
    corrected, meta = posse.align(ref, target)
    
    print(f"\nPOSSEv6 correction:")
    print(f"  Alpha (effective): {meta['alpha_mean']:.4f}")
    print(f"  Beta (effective): {meta['beta_mean']:.4f}")
    print(f"  Alpha (internal): {meta['alpha_internal']:.4f}")
    print(f"  Beta (internal): {meta['beta_internal']:.4f}")
    
    print(f"\nCorrected data:")
    print(f"  Mean: {np.mean(corrected.data):.3f}, Std: {np.std(corrected.data):.3f}")
    print(f"  Range: [{np.min(corrected.data):.3f}, {np.max(corrected.data):.3f}]")
    
    # What GMM would do (simple moment matching)
    ref_mean = np.mean(train_data, axis=1, keepdims=True)
    ref_std = np.std(train_data, axis=1, keepdims=True) + 1e-8
    target_mean = np.mean(test_data, axis=1, keepdims=True)
    target_std = np.std(test_data, axis=1, keepdims=True) + 1e-8
    
    gmm_alpha = ref_std / target_std
    gmm_beta = ref_mean - gmm_alpha * target_mean
    gmm_corrected = gmm_alpha * test_data + gmm_beta
    
    print(f"\nGMM (moment matching) correction:")
    print(f"  Alpha mean: {np.mean(gmm_alpha):.4f}")
    print(f"  Beta mean: {np.mean(gmm_beta):.4f}")
    print(f"  Corrected mean: {np.mean(gmm_corrected):.3f}, Std: {np.std(gmm_corrected):.3f}")
    
    # The problem: POSSEv6 isn't correcting enough!
    print(f"\n" + "=" * 60)
    print("DIAGNOSIS:")
    print("=" * 60)
    
    # Check if corrected data matches reference distribution
    ref_global_mean = np.mean(train_data)
    ref_global_std = np.std(train_data)
    
    posse_mean_error = abs(np.mean(corrected.data) - ref_global_mean)
    posse_std_error = abs(np.std(corrected.data) - ref_global_std)
    
    gmm_mean_error = abs(np.mean(gmm_corrected) - ref_global_mean)
    gmm_std_error = abs(np.std(gmm_corrected) - ref_global_std)
    
    print(f"Reference: mean={ref_global_mean:.3f}, std={ref_global_std:.3f}")
    print(f"\nPOSSEv6 alignment error:")
    print(f"  Mean error: {posse_mean_error:.3f}")
    print(f"  Std error: {posse_std_error:.3f}")
    
    print(f"\nGMM alignment error:")
    print(f"  Mean error: {gmm_mean_error:.3f}")
    print(f"  Std error: {gmm_std_error:.3f}")

if __name__ == "__main__":
    main()
