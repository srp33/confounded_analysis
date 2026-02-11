#!/usr/bin/env python3
"""Cross-study prediction test - the real challenge"""

import sys
sys.path.append('scripts')

import numpy as np
import subprocess
import os
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, roc_auc_score

from pace import BatchData
from posse import POSSE, POSSEHyperparameters
from posse_v6 import POSSEv6, POSSEv6Hyperparameters, BatchData as BatchDataV6

def load_two_studies():
    """Load two studies from R"""
    r_script = """
    load('data/TB_real_data.RData')
    
    # Study 1: Africa (training)
    write.table(dat_lst[['Africa']], 'temp_train_data.txt', sep='\\t', row.names=FALSE, col.names=FALSE)
    write.table(label_lst[['Africa']], 'temp_train_labels.txt', row.names=FALSE, col.names=FALSE)
    
    # Study 2: GSE37250_M (test)
    write.table(dat_lst[['GSE37250_M']], 'temp_test_data.txt', sep='\\t', row.names=FALSE, col.names=FALSE)
    write.table(label_lst[['GSE37250_M']], 'temp_test_labels.txt', row.names=FALSE, col.names=FALSE)
    
    # Gene names
    write.table(rownames(dat_lst[['Africa']]), 'temp_genes.txt', row.names=FALSE, col.names=FALSE, quote=FALSE)
    """
    subprocess.run(['R', '--slave', '--vanilla'], input=r_script, text=True, capture_output=True)
    
    train_data = np.loadtxt('temp_train_data.txt', delimiter='\t')
    train_labels = np.loadtxt('temp_train_labels.txt')
    test_data = np.loadtxt('temp_test_data.txt', delimiter='\t')
    test_labels = np.loadtxt('temp_test_labels.txt')
    genes = np.loadtxt('temp_genes.txt', dtype=str)
    
    for f in ['temp_train_data.txt', 'temp_train_labels.txt', 'temp_test_data.txt', 'temp_test_labels.txt', 'temp_genes.txt']:
        if os.path.exists(f):
            os.remove(f)
    
    return train_data, train_labels, test_data, test_labels, genes

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

def run_naive(ref_data, target_data):
    """Naive correction: global moment matching (single alpha/beta for all genes)"""
    ref_mean = np.mean(ref_data)
    ref_std = np.std(ref_data) + 1e-8
    target_mean = np.mean(target_data)
    target_std = np.std(target_data) + 1e-8
    
    alpha = ref_std / target_std
    beta = ref_mean - alpha * target_mean
    return alpha * target_data + beta

def run_gmm(ref_data, target_data):
    """GMM-based correction: fit mixture models, then gene-wise moment matching"""
    from sklearn.mixture import GaussianMixture
    
    # Fit GMM to identify clusters in each dataset
    ref_samples = ref_data.T  # samples x genes
    target_samples = target_data.T
    
    # Use subset of genes for GMM fitting
    n_genes_subset = min(500, ref_samples.shape[1])
    np.random.seed(42)
    gene_idx = np.random.choice(ref_samples.shape[1], n_genes_subset, replace=False)
    
    # Fit GMMs to identify sample clusters
    gmm_ref = GaussianMixture(n_components=2, random_state=42, n_init=3)
    gmm_target = GaussianMixture(n_components=2, random_state=42, n_init=3)
    
    gmm_ref.fit(ref_samples[:, gene_idx])
    gmm_target.fit(target_samples[:, gene_idx])
    
    # Get cluster assignments
    ref_clusters = gmm_ref.predict(ref_samples[:, gene_idx])
    target_clusters = gmm_target.predict(target_samples[:, gene_idx])
    
    # Compute cluster-aware gene-wise correction
    # For each gene, match moments considering cluster structure
    corrected = np.zeros_like(target_data)
    
    for g in range(ref_data.shape[0]):
        # Weighted mean/std based on cluster proportions
        ref_gene = ref_data[g, :]
        target_gene = target_data[g, :]
        
        # Cluster-weighted statistics for reference
        ref_c0_mask = ref_clusters == 0
        ref_c1_mask = ref_clusters == 1
        
        if np.sum(ref_c0_mask) > 1 and np.sum(ref_c1_mask) > 1:
            ref_mean_c0 = np.mean(ref_gene[ref_c0_mask])
            ref_mean_c1 = np.mean(ref_gene[ref_c1_mask])
            ref_std_c0 = np.std(ref_gene[ref_c0_mask]) + 1e-8
            ref_std_c1 = np.std(ref_gene[ref_c1_mask]) + 1e-8
            
            # Target cluster stats
            target_c0_mask = target_clusters == 0
            target_c1_mask = target_clusters == 1
            
            if np.sum(target_c0_mask) > 1 and np.sum(target_c1_mask) > 1:
                target_mean_c0 = np.mean(target_gene[target_c0_mask])
                target_mean_c1 = np.mean(target_gene[target_c1_mask])
                target_std_c0 = np.std(target_gene[target_c0_mask]) + 1e-8
                target_std_c1 = np.std(target_gene[target_c1_mask]) + 1e-8
                
                # Correct each cluster separately
                corrected_gene = target_gene.copy()
                
                # Cluster 0
                alpha_c0 = ref_std_c0 / target_std_c0
                beta_c0 = ref_mean_c0 - alpha_c0 * target_mean_c0
                corrected_gene[target_c0_mask] = alpha_c0 * target_gene[target_c0_mask] + beta_c0
                
                # Cluster 1
                alpha_c1 = ref_std_c1 / target_std_c1
                beta_c1 = ref_mean_c1 - alpha_c1 * target_mean_c1
                corrected_gene[target_c1_mask] = alpha_c1 * target_gene[target_c1_mask] + beta_c1
                
                corrected[g, :] = corrected_gene
                continue
        
        # Fallback to simple gene-wise moment matching
        ref_mean_g = np.mean(ref_gene)
        ref_std_g = np.std(ref_gene) + 1e-8
        target_mean_g = np.mean(target_gene)
        target_std_g = np.std(target_gene) + 1e-8
        
        alpha = ref_std_g / target_std_g
        beta = ref_mean_g - alpha * target_mean_g
        corrected[g, :] = alpha * target_gene + beta
    
    return corrected

def evaluate_classifier(train_data, train_labels, test_data, test_labels, method_name):
    """Train SVM on train, evaluate on test"""
    # Transpose: genes x samples -> samples x genes
    X_train = train_data.T
    X_test = test_data.T
    
    # Standardize
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train SVM
    clf = SVC(kernel='linear', probability=True, random_state=42)
    clf.fit(X_train_scaled, train_labels)
    
    # Predict
    y_pred = clf.predict(X_test_scaled)
    y_prob = clf.predict_proba(X_test_scaled)[:, 1]
    
    acc = accuracy_score(test_labels, y_pred)
    try:
        auc = roc_auc_score(test_labels, y_prob)
    except:
        auc = 0.5
    
    return acc, auc

def main():
    print("CROSS-STUDY TB PREDICTION TEST")
    print("=" * 60)
    print("Train: Africa study | Test: GSE37250_M study")
    print("=" * 60)
    
    train_data, train_labels, test_data, test_labels, genes = load_two_studies()
    pathways = get_pathways(genes)
    
    print(f"Train: {train_data.shape[1]} samples, Test: {test_data.shape[1]} samples")
    print(f"Genes: {train_data.shape[0]}")
    print()
    
    results = {}
    
    # No correction (baseline)
    print("1. NO CORRECTION (baseline)")
    acc, auc = evaluate_classifier(train_data, train_labels, test_data, test_labels, "None")
    results['None'] = {'acc': acc, 'auc': auc}
    print(f"   Accuracy: {acc:.3f}, AUC: {auc:.3f}")
    
    # Naive (global moment matching)
    print("\n2. NAIVE (global moment matching)")
    try:
        corrected = run_naive(train_data, test_data)
        acc, auc = evaluate_classifier(train_data, train_labels, corrected, test_labels, "Naive")
        results['Naive'] = {'acc': acc, 'auc': auc}
        print(f"   Accuracy: {acc:.3f}, AUC: {auc:.3f}")
    except Exception as e:
        print(f"   CRASHED: {e}")
        results['Naive'] = {'acc': 0, 'auc': 0.5}
    
    # POSSEv5
    print("\n3. POSSEv5")
    try:
        ref = BatchData(data=train_data, gene_indices=genes)
        target = BatchData(data=test_data, gene_indices=genes)
        posse = POSSE(pathway_dict=pathways, hyperparams=POSSEHyperparameters(tau=25.0))
        corrected, meta = posse.align(ref, target)
        acc, auc = evaluate_classifier(train_data, train_labels, corrected.data, test_labels, "POSSEv5")
        results['POSSEv5'] = {'acc': acc, 'auc': auc}
        print(f"   Accuracy: {acc:.3f}, AUC: {auc:.3f}")
    except Exception as e:
        print(f"   CRASHED: {e}")
        results['POSSEv5'] = {'acc': 0, 'auc': 0.5}
    
    # POSSEv6
    print("\n4. POSSEv6")
    try:
        ref = BatchDataV6(data=train_data, gene_indices=genes)
        target = BatchDataV6(data=test_data, gene_indices=genes)
        posse = POSSEv6(pathway_dict=pathways, hyperparams=POSSEv6Hyperparameters(tau=25.0))
        corrected, meta = posse.align(ref, target)
        acc, auc = evaluate_classifier(train_data, train_labels, corrected.data, test_labels, "POSSEv6")
        results['POSSEv6'] = {'acc': acc, 'auc': auc}
        print(f"   Accuracy: {acc:.3f}, AUC: {auc:.3f}")
    except Exception as e:
        print(f"   CRASHED: {e}")
        results['POSSEv6'] = {'acc': 0, 'auc': 0.5}
    
    # ComBat
    print("\n5. ComBat")
    try:
        corrected = run_combat(train_data, test_data)
        acc, auc = evaluate_classifier(train_data, train_labels, corrected, test_labels, "ComBat")
        results['ComBat'] = {'acc': acc, 'auc': auc}
        print(f"   Accuracy: {acc:.3f}, AUC: {auc:.3f}")
    except Exception as e:
        print(f"   CRASHED: {e}")
        results['ComBat'] = {'acc': 0, 'auc': 0.5}
    
    # GMM (cluster-aware correction)
    print("\n6. GMM (cluster-aware correction)")
    try:
        corrected = run_gmm(train_data, test_data)
        acc, auc = evaluate_classifier(train_data, train_labels, corrected, test_labels, "GMM")
        results['GMM'] = {'acc': acc, 'auc': auc}
        print(f"   Accuracy: {acc:.3f}, AUC: {auc:.3f}")
    except Exception as e:
        print(f"   CRASHED: {e}")
        results['GMM'] = {'acc': 0, 'auc': 0.5}
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY (ranked by AUC)")
    print("=" * 60)
    sorted_results = sorted(results.items(), key=lambda x: x[1]['auc'], reverse=True)
    for i, (method, res) in enumerate(sorted_results, 1):
        print(f"  {i}. {method}: AUC={res['auc']:.3f}, Acc={res['acc']:.3f}")
    
    # Improvement over baseline
    baseline_auc = results['None']['auc']
    print(f"\nImprovement over baseline (AUC={baseline_auc:.3f}):")
    for method, res in sorted_results:
        if method != 'None':
            delta = res['auc'] - baseline_auc
            sign = '+' if delta >= 0 else ''
            print(f"  {method}: {sign}{delta:.3f}")

if __name__ == "__main__":
    main()
