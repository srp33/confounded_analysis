#!/usr/bin/env python3
"""Quick comparative test - POSSEv6 vs ComBat vs GMM"""

import sys
sys.path.append('scripts')

import numpy as np
import subprocess
import os

from pace import BatchData
from posse import POSSE, POSSEHyperparameters
from posse_v6 import POSSEv6, POSSEv6Hyperparameters, BatchData as BatchDataV6

def load_one_study():
    """Load one study from R"""
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

def run_combat(ref_data, target_data):
    """Run ComBat via R"""
    np.savetxt('temp_ref.txt', ref_data, delimiter='\t')
    np.savetxt('temp_target.txt', target_data, delimiter='\t')
    
    r_script = """
    library(sva)
    ref <- read.table('temp_ref.txt', sep='\\t')
    target <- read.table('temp_target.txt', sep='\\t')
    combined <- cbind(ref, target)
    batch <- c(rep(1, ncol(ref)), rep(2, ncol(target)))
    sink('/dev/null')
    corrected <- ComBat(dat=as.matrix(combined), batch=batch)
    sink()
    corrected_target <- corrected[, (ncol(ref)+1):ncol(corrected)]
    write.table(corrected_target, 'temp_combat_out.txt', sep='\\t', row.names=FALSE, col.names=FALSE)
    """
    result = subprocess.run(['R', '--slave', '--vanilla'], input=r_script, text=True, capture_output=True)
    
    if result.returncode != 0:
        raise Exception(f"ComBat failed: {result.stderr}")
    
    corrected = np.loadtxt('temp_combat_out.txt', delimiter='\t')
    for f in ['temp_ref.txt', 'temp_target.txt', 'temp_combat_out.txt']:
        if os.path.exists(f):
            os.remove(f)
    return corrected

def run_gmm(ref_data, target_data):
    """Simple GMM-based correction"""
    ref_mean = np.mean(ref_data, axis=1, keepdims=True)
    ref_std = np.std(ref_data, axis=1, keepdims=True) + 1e-8
    target_mean = np.mean(target_data, axis=1, keepdims=True)
    target_std = np.std(target_data, axis=1, keepdims=True) + 1e-8
    
    alpha = ref_std / target_std
    beta = ref_mean - alpha * target_mean
    return alpha * target_data + beta

def test_silence():
    """Test: identical data should return identical data"""
    print("=" * 60)
    print("TEST: SILENCE (identical data → α=1, β=0)")
    print("=" * 60)
    
    data, genes = load_one_study()
    pathways = get_pathways(genes)
    
    results = {}
    
    # POSSEv5
    try:
        ref = BatchData(data=data, gene_indices=genes)
        target = BatchData(data=data.copy(), gene_indices=genes)
        posse = POSSE(pathway_dict=pathways, hyperparams=POSSEHyperparameters(tau=25.0))
        corrected, meta = posse.align(ref, target)
        rmse = np.sqrt(np.mean((corrected.data - data)**2))
        results['POSSEv5'] = {'rmse': rmse, 'alpha': meta.get('alpha_final_mean', 1.0)}
        print(f"  POSSEv5:  RMSE={rmse:.6f}, α={results['POSSEv5']['alpha']:.4f}")
    except Exception as e:
        print(f"  POSSEv5:  CRASHED - {e}")
        results['POSSEv5'] = {'rmse': float('inf'), 'alpha': None}
    
    # POSSEv6
    try:
        ref = BatchDataV6(data=data, gene_indices=genes)
        target = BatchDataV6(data=data.copy(), gene_indices=genes)
        posse = POSSEv6(pathway_dict=pathways, hyperparams=POSSEv6Hyperparameters(tau=25.0))
        corrected, meta = posse.align(ref, target)
        rmse = np.sqrt(np.mean((corrected.data - data)**2))
        results['POSSEv6'] = {'rmse': rmse, 'alpha': meta.get('alpha_mean', 1.0)}
        print(f"  POSSEv6:  RMSE={rmse:.6f}, α={results['POSSEv6']['alpha']:.4f}")
    except Exception as e:
        print(f"  POSSEv6:  CRASHED - {e}")
        results['POSSEv6'] = {'rmse': float('inf'), 'alpha': None}
    
    # ComBat
    try:
        corrected = run_combat(data, data.copy())
        rmse = np.sqrt(np.mean((corrected - data)**2))
        results['ComBat'] = {'rmse': rmse, 'alpha': 1.0}
        print(f"  ComBat:   RMSE={rmse:.6f}")
    except Exception as e:
        print(f"  ComBat:   CRASHED - {e}")
        results['ComBat'] = {'rmse': float('inf'), 'alpha': None}
    
    # GMM
    try:
        corrected = run_gmm(data, data.copy())
        rmse = np.sqrt(np.mean((corrected - data)**2))
        results['GMM'] = {'rmse': rmse, 'alpha': 1.0}
        print(f"  GMM:      RMSE={rmse:.6f}")
    except Exception as e:
        print(f"  GMM:      CRASHED - {e}")
        results['GMM'] = {'rmse': float('inf'), 'alpha': None}
    
    return results

def test_gene_specific():
    """Test: gene-specific batch effects (realistic)"""
    print("\n" + "=" * 60)
    print("TEST: GENE-SPECIFIC BATCH EFFECTS (realistic)")
    print("=" * 60)
    
    data, genes = load_one_study()
    pathways = get_pathways(genes)
    
    # Create gene-specific batch effects
    n_genes = data.shape[0]
    np.random.seed(42)
    gene_alpha = np.random.normal(1.2, 0.2, size=(n_genes, 1))
    gene_beta = np.random.normal(0.5, 0.5, size=(n_genes, 1))
    distorted = data * gene_alpha + gene_beta
    
    print(f"  Applied: α~N(1.2,0.2), β~N(0.5,0.5) per gene")
    
    results = {}
    
    # POSSEv5
    try:
        ref = BatchData(data=data, gene_indices=genes)
        target = BatchData(data=distorted, gene_indices=genes)
        posse = POSSE(pathway_dict=pathways, hyperparams=POSSEHyperparameters(tau=25.0))
        corrected, meta = posse.align(ref, target)
        rmse = np.sqrt(np.mean((corrected.data - data)**2))
        results['POSSEv5'] = {'rmse': rmse}
        print(f"  POSSEv5:  Recovery RMSE={rmse:.4f}")
    except Exception as e:
        print(f"  POSSEv5:  CRASHED - {e}")
        results['POSSEv5'] = {'rmse': float('inf')}
    
    # POSSEv6
    try:
        ref = BatchDataV6(data=data, gene_indices=genes)
        target = BatchDataV6(data=distorted, gene_indices=genes)
        posse = POSSEv6(pathway_dict=pathways, hyperparams=POSSEv6Hyperparameters(tau=25.0))
        corrected, meta = posse.align(ref, target)
        rmse = np.sqrt(np.mean((corrected.data - data)**2))
        results['POSSEv6'] = {'rmse': rmse}
        print(f"  POSSEv6:  Recovery RMSE={rmse:.4f}")
    except Exception as e:
        print(f"  POSSEv6:  CRASHED - {e}")
        results['POSSEv6'] = {'rmse': float('inf')}
    
    # ComBat
    try:
        corrected = run_combat(data, distorted)
        rmse = np.sqrt(np.mean((corrected - data)**2))
        results['ComBat'] = {'rmse': rmse}
        print(f"  ComBat:   Recovery RMSE={rmse:.4f}")
    except Exception as e:
        print(f"  ComBat:   CRASHED - {e}")
        results['ComBat'] = {'rmse': float('inf')}
    
    # GMM
    try:
        corrected = run_gmm(data, distorted)
        rmse = np.sqrt(np.mean((corrected - data)**2))
        results['GMM'] = {'rmse': rmse}
        print(f"  GMM:      Recovery RMSE={rmse:.4f}")
    except Exception as e:
        print(f"  GMM:      CRASHED - {e}")
        results['GMM'] = {'rmse': float('inf')}
    
    return results

if __name__ == "__main__":
    print("QUICK COMPARATIVE TEST: POSSEv5 vs POSSEv6 vs ComBat vs GMM")
    print("=" * 60)
    
    silence_results = test_silence()
    gene_specific_results = test_gene_specific()
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    # Rank by silence test
    print("\nSilence Test (lower RMSE = better):")
    sorted_silence = sorted(silence_results.items(), key=lambda x: x[1]['rmse'])
    for i, (method, res) in enumerate(sorted_silence, 1):
        print(f"  {i}. {method}: RMSE={res['rmse']:.6f}")
    
    # Rank by gene-specific test
    print("\nGene-Specific Test (lower RMSE = better):")
    sorted_gene = sorted(gene_specific_results.items(), key=lambda x: x[1]['rmse'])
    for i, (method, res) in enumerate(sorted_gene, 1):
        print(f"  {i}. {method}: RMSE={res['rmse']:.4f}")
