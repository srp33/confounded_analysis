#!/usr/bin/env python3
"""Focused POSSE debugging - trace exactly where it fails"""

import sys
sys.path.append('scripts')

import numpy as np
import subprocess
import os

from pace import BatchData
from posse import POSSE, POSSEHyperparameters

class SuppressOutput:
    def __init__(self):
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for _ in range(2)]
        self.save_fds = [os.dup(1), os.dup(2)]
    def __enter__(self):
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)
    def __exit__(self, *_):
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        for fd in self.null_fds + self.save_fds:
            os.close(fd)

def load_one_study():
    """Load just one study for debugging"""
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
    """Get minimal pathway dict"""
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
    """Test 1: Identical data should give α=1, β=0"""
    print("=" * 60)
    print("TEST 1: SILENCE TEST (identical data)")
    print("=" * 60)
    
    data, genes = load_one_study()
    pathways = get_pathways(genes)
    
    ref = BatchData(data=data, gene_indices=genes)
    target = BatchData(data=data.copy(), gene_indices=genes)
    
    posse = POSSE(pathway_dict=pathways, hyperparams=POSSEHyperparameters(tau=25.0, top_k_percent=0.2))
    
    with SuppressOutput():
        corrected, meta = posse.align(ref, target)
    
    rmse = np.sqrt(np.mean((corrected.data - data)**2))
    print(f"RMSE (should be ~0): {rmse:.6f}")
    print(f"Alpha mean (should be 1.0): {meta.get('alpha_mean', 'N/A')}")
    print(f"Beta mean (should be 0.0): {meta.get('beta_mean', 'N/A')}")
    
    # Check what's happening internally
    print("\nDiagnostic history:")
    for h in meta.get('diagnostic_history', [])[:2]:
        print(f"  Iter {h['iter']}: α_range={h['alpha_range']}, β_range={h['beta_range']}")
    
    return rmse < 0.01

def test_known_shift():
    """Test 2: Known shift should be recovered"""
    print("\n" + "=" * 60)
    print("TEST 2: KNOWN SHIFT (Y = X + 2.0)")
    print("=" * 60)
    
    data, genes = load_one_study()
    pathways = get_pathways(genes)
    
    # Apply known shift
    shift = 2.0
    shifted_data = data + shift
    
    ref = BatchData(data=data, gene_indices=genes)
    target = BatchData(data=shifted_data, gene_indices=genes)
    
    posse = POSSE(pathway_dict=pathways, hyperparams=POSSEHyperparameters(tau=25.0, top_k_percent=0.2))
    
    with SuppressOutput():
        corrected, meta = posse.align(ref, target)
    
    # POSSE should apply β ≈ -2.0 to recover original
    recovery_rmse = np.sqrt(np.mean((corrected.data - data)**2))
    print(f"Recovery RMSE (should be ~0): {recovery_rmse:.6f}")
    print(f"Alpha mean (should be ~1.0): {meta.get('alpha_mean', 'N/A')}")
    print(f"Beta mean (should be ~-2.0): {meta.get('beta_mean', 'N/A')}")
    
    # What did POSSE actually estimate?
    print("\nDiagnostic history:")
    for h in meta.get('diagnostic_history', [])[:2]:
        print(f"  Iter {h['iter']}: α_range={h['alpha_range']}, β_range={h['beta_range']}")
    
    return recovery_rmse < 0.5

def test_known_scale():
    """Test 3: Known scale should be recovered"""
    print("\n" + "=" * 60)
    print("TEST 3: KNOWN SCALE (Y = 2.0 * X)")
    print("=" * 60)
    
    data, genes = load_one_study()
    pathways = get_pathways(genes)
    
    # Apply known scale
    scale = 2.0
    scaled_data = data * scale
    
    ref = BatchData(data=data, gene_indices=genes)
    target = BatchData(data=scaled_data, gene_indices=genes)
    
    posse = POSSE(pathway_dict=pathways, hyperparams=POSSEHyperparameters(tau=25.0, top_k_percent=0.2))
    
    with SuppressOutput():
        corrected, meta = posse.align(ref, target)
    
    # POSSE should apply α ≈ 0.5 to recover original
    recovery_rmse = np.sqrt(np.mean((corrected.data - data)**2))
    print(f"Recovery RMSE (should be ~0): {recovery_rmse:.6f}")
    print(f"Alpha mean (should be ~0.5): {meta.get('alpha_mean', 'N/A')}")
    print(f"Beta mean (should be ~0.0): {meta.get('beta_mean', 'N/A')}")
    
    return recovery_rmse < 0.5

def trace_combat_priors():
    """Trace what ComBat priors are being computed"""
    print("\n" + "=" * 60)
    print("TRACE: COMBAT PRIOR COMPUTATION")
    print("=" * 60)
    
    data, genes = load_one_study()
    
    # Create identical data scenario
    X = data
    Y = data.copy()
    
    # Manually trace ComBat computation
    from posse import ComBatBaseline
    combat = ComBatBaseline()
    
    # Suppress the verbose output but capture the priors
    with SuppressOutput():
        alpha_prior, beta_prior = combat.compute_baseline(X, Y)
    
    print(f"For IDENTICAL data:")
    print(f"  Alpha prior: mean={np.mean(alpha_prior):.4f}, std={np.std(alpha_prior):.4f}")
    print(f"  Alpha range: [{np.min(alpha_prior):.4f}, {np.max(alpha_prior):.4f}]")
    print(f"  Beta prior: mean={np.mean(beta_prior):.4f}, std={np.std(beta_prior):.4f}")
    print(f"  Beta range: [{np.min(beta_prior):.4f}, {np.max(beta_prior):.4f}]")
    
    # Now with shifted data
    Y_shifted = data + 2.0
    with SuppressOutput():
        alpha_prior_s, beta_prior_s = combat.compute_baseline(X, Y_shifted)
    
    print(f"\nFor SHIFTED data (Y = X + 2.0):")
    print(f"  Alpha prior: mean={np.mean(alpha_prior_s):.4f}")
    print(f"  Beta prior: mean={np.mean(beta_prior_s):.4f} (should be ~-2.0)")

def trace_preprocessing():
    """Trace what preprocessing does to the data"""
    print("\n" + "=" * 60)
    print("TRACE: PREPROCESSING EFFECT")
    print("=" * 60)
    
    data, genes = load_one_study()
    pathways = get_pathways(genes)
    
    posse = POSSE(pathway_dict=pathways, hyperparams=POSSEHyperparameters())
    
    X = data
    Y = data + 2.0  # Known shift
    
    print(f"Before preprocessing:")
    print(f"  X mean: {np.mean(X):.4f}, Y mean: {np.mean(Y):.4f}")
    print(f"  Difference: {np.mean(Y) - np.mean(X):.4f} (should be 2.0)")
    
    # Call preprocessing
    X_p, Y_p = posse.adaptive_preprocessing(X, Y)
    
    print(f"\nAfter preprocessing (asinh transform):")
    print(f"  X_p mean: {np.mean(X_p):.4f}, Y_p mean: {np.mean(Y_p):.4f}")
    print(f"  Difference: {np.mean(Y_p) - np.mean(X_p):.4f}")
    
    # Check if preprocessing destroys the shift information
    print(f"\nPreprocessing effect on shift:")
    print(f"  Original shift: 2.0")
    print(f"  Preprocessed shift: {np.mean(Y_p) - np.mean(X_p):.4f}")

if __name__ == "__main__":
    print("POSSE DEBUGGING SESSION")
    print("=" * 60)
    
    trace_preprocessing()
    trace_combat_priors()
    
    print("\n" + "=" * 60)
    print("RUNNING VALIDATION TESTS")
    print("=" * 60)
    
    silence_pass = test_silence()
    shift_pass = test_known_shift()
    scale_pass = test_known_scale()
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Silence test: {'PASS' if silence_pass else 'FAIL'}")
    print(f"Shift recovery: {'PASS' if shift_pass else 'FAIL'}")
    print(f"Scale recovery: {'PASS' if scale_pass else 'FAIL'}")
