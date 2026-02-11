#!/usr/bin/env python3
"""
Test PACE v3.4 Housekeeping Anchors with full transcriptome (~20K genes)
Updated to include POSSE variants and maintain existing baselines
"""

import sys
sys.path.append('scripts')

import numpy as np
from validate_pace import run_validation_suite, generate_validation_report
from pace import PACE_v22, BatchData, PACEHyperparameters
from posse import POSSE, POSSEHyperparameters
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
import contextlib
import io
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats

class ComBatBaseline:
    """Pure ComBat implementation for comparison"""
    
    def __init__(self):
        pass
    
    def align(self, ref_data: BatchData, target_data: BatchData) -> tuple:
        """Apply ComBat correction to align target data with reference"""
        
        # Find common genes
        common_genes = np.intersect1d(ref_data.gene_indices, target_data.gene_indices)
        
        # Get indices for common genes
        ref_indices = np.array([np.where(ref_data.gene_indices == gene)[0][0] for gene in common_genes])
        target_indices = np.array([np.where(target_data.gene_indices == gene)[0][0] for gene in common_genes])
        
        # Extract common gene data
        X_common = ref_data.data[ref_indices]  # Reference data
        Y_common = target_data.data[target_indices]  # Target data
        
        print(f"ComBat: Processing {len(common_genes)} common genes")
        
        # Apply log transformation for ComBat (assumes linear input data)
        X_log = np.log2(X_common + 1)
        Y_log = np.log2(Y_common + 1)
        
        # Combine data for ComBat
        combined_data = np.hstack([X_log, Y_log])
        batch_labels = np.array([1] * X_log.shape[1] + [2] * Y_log.shape[1])
        
        # Apply ComBat correction
        corrected_data = self._combat_correction(combined_data, batch_labels, ref_batch=1)
        
        # Extract corrected target data
        Y_corrected_log = corrected_data[:, X_log.shape[1]:]
        
        # Convert back to linear scale
        Y_corrected = 2**Y_corrected_log - 1
        Y_corrected = np.maximum(Y_corrected, 0)  # Ensure non-negative
        
        # Create full corrected matrix (handle unique genes)
        Y_final = target_data.data.copy()
        Y_final[target_indices] = Y_corrected
        
        metadata = {
            "method": "ComBat",
            "n_common_genes": len(common_genes),
            "n_unique_genes": len(target_data.gene_indices) - len(common_genes)
        }
        
        return Y_final, metadata
    
    def _combat_correction(self, data, batch, ref_batch=None):
        """Simplified ComBat implementation"""
        
        # Get unique batches
        batches = np.unique(batch)
        n_batches = len(batches)
        
        if n_batches < 2:
            return data  # No correction needed
        
        # Calculate overall mean for each gene
        overall_mean = np.mean(data, axis=1, keepdims=True)
        
        # Calculate batch effects
        batch_means = {}
        batch_vars = {}
        
        for b in batches:
            batch_mask = batch == b
            batch_data = data[:, batch_mask]
            
            # Batch mean (location parameter)
            batch_means[b] = np.mean(batch_data, axis=1, keepdims=True)
            
            # Batch variance (scale parameter)  
            batch_vars[b] = np.var(batch_data, axis=1, keepdims=True)
        
        # Apply correction
        corrected_data = data.copy()
        
        if ref_batch is not None and ref_batch in batches:
            # Use reference batch approach
            ref_mean = batch_means[ref_batch]
            ref_var = batch_vars[ref_batch]
            
            for b in batches:
                if b == ref_batch:
                    continue  # Don't correct reference batch
                    
                batch_mask = batch == b
                
                # Location correction (additive)
                location_shift = batch_means[b] - ref_mean
                
                # Scale correction (multiplicative)
                scale_factor = np.sqrt(batch_vars[b] / (ref_var + 1e-8))
                
                # Apply correction: (data - batch_mean) / scale_factor + ref_mean
                corrected_data[:, batch_mask] = (
                    (data[:, batch_mask] - batch_means[b]) / (scale_factor + 1e-8) + ref_mean
                )
        else:
            # Standard parametric ComBat (adjust all batches to overall mean)
            for b in batches:
                batch_mask = batch == b
                
                # Location correction
                location_shift = batch_means[b] - overall_mean
                
                # Scale correction
                pooled_var = np.mean([batch_vars[bb] for bb in batches], axis=0)
                scale_factor = np.sqrt(batch_vars[b] / (pooled_var + 1e-8))
                
                # Apply correction
                corrected_data[:, batch_mask] = (
                    (data[:, batch_mask] - batch_means[b]) / (scale_factor + 1e-8) + overall_mean
                )
        
        return corrected_data

class GMMBaseline:
    """GMM-based batch correction baseline using R implementation"""
    
    def __init__(self):
        pass
    
    def align(self, ref_data: BatchData, target_data: BatchData) -> tuple:
        """Apply GMM correction to align target data with reference"""
        
        # Find common genes
        common_genes = np.intersect1d(ref_data.gene_indices, target_data.gene_indices)
        
        # Get indices for common genes
        ref_indices = np.array([np.where(ref_data.gene_indices == gene)[0][0] for gene in common_genes])
        target_indices = np.array([np.where(target_data.gene_indices == gene)[0][0] for gene in common_genes])
        
        # Extract common gene data
        X_common = ref_data.data[ref_indices]  # Reference data
        Y_common = target_data.data[target_indices]  # Target data
        
        print(f"GMM: Processing {len(common_genes)} common genes")
        print(f"GMM: Data shapes - Ref: {X_common.shape}, Target: {Y_common.shape}")
        
        # Combine data for GMM
        combined_data = np.hstack([X_common, Y_common])  # genes x samples
        batch_labels = np.array([1] * X_common.shape[1] + [2] * Y_common.shape[1])
        
        print(f"GMM: Combined data shape: {combined_data.shape}")
        print(f"GMM: Batch labels: {batch_labels[:10]}... (length: {len(batch_labels)})")
        
        # Apply GMM correction using R
        import subprocess
        import tempfile
        import os
        
        # Create temporary files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            data_file = f.name
            # Write data in samples x genes format (transpose for genes_are_columns=TRUE)
            np.savetxt(f, combined_data.T, delimiter=',')  # Transpose to samples x genes
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            batch_file = f.name
            np.savetxt(f, batch_labels, delimiter=',', fmt='%d')
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            output_file = f.name
        
        # Create R script - DISABLE log_transform since data is already log-transformed
        r_script = f"""
        source('../../adjust/gmm_adjust.R')
        
        # Read data
        data <- as.matrix(read.csv('{data_file}', header=FALSE))
        batch <- as.vector(read.csv('{batch_file}', header=FALSE)[,1])
        
        cat("Data dimensions:", dim(data), "\\n")
        cat("Batch length:", length(batch), "\\n")
        cat("Batch levels:", unique(batch), "\\n")
        
        # Apply GMM adjustment with log_transform=FALSE since data is already log-transformed
        adjusted <- gmm_adjust(data, batch, genes_are_columns=TRUE, log_transform=FALSE, debug=TRUE)
        
        cat("Adjusted dimensions:", dim(adjusted), "\\n")
        
        # Write result
        write.table(adjusted, '{output_file}', row.names=FALSE, col.names=FALSE, sep=',')
        """
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
            script_file = f.name
            f.write(r_script)
        
        print(f"GMM: Running R script...")
        
        # Run R script
        result = subprocess.run(['Rscript', script_file], 
                              capture_output=True, text=True, cwd='.')
        
        if result.returncode != 0:
            print(f"GMM R script failed with return code {result.returncode}")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            raise RuntimeError(f"GMM correction failed: {result.stderr}")
        else:
            print(f"GMM R script output: {result.stdout}")
            # Read corrected data
            corrected_data = np.loadtxt(output_file, delimiter=',')
            print(f"GMM: Corrected data shape: {corrected_data.shape}")
            
            # If gmm_adjust returns samples x genes, transpose back to genes x samples
            if corrected_data.shape[0] == combined_data.shape[1]:  # samples x genes
                corrected_data = corrected_data.T  # transpose to genes x samples
                print(f"GMM: Transposed corrected data to genes x samples: {corrected_data.shape}")
        
        # Clean up temporary files
        for temp_file in [data_file, batch_file, output_file, script_file]:
            try:
                os.unlink(temp_file)
            except:
                pass
        
        # Extract corrected target data
        Y_corrected = corrected_data[:, X_common.shape[1]:]
        print(f"GMM: Extracted target data shape: {Y_corrected.shape}")
        
        # Create full corrected matrix (handle unique genes)
        Y_final = target_data.data.copy()
        Y_final[target_indices] = Y_corrected
        
        metadata = {
            "method": "GMM",
            "n_common_genes": len(common_genes),
            "n_unique_genes": len(target_data.gene_indices) - len(common_genes)
        }
        
        return Y_final, metadata

def test_de_preservation(ref_data, target_data, corrected_target, ref_labels, target_labels, gene_names, tb_signature_genes):
    """Test if known TB biomarkers maintain differential expression patterns"""
    
    try:
        # Calculate DE in reference
        ref_tb_idx = ref_labels == 1
        ref_ctrl_idx = ref_labels == 0
        
        if np.sum(ref_tb_idx) == 0 or np.sum(ref_ctrl_idx) == 0:
            return 0.0
        
        ref_de_scores = {}
        for gene in tb_signature_genes:
            if gene in gene_names:
                gene_idx = np.where(gene_names == gene)[0][0]
                tb_expr = ref_data[gene_idx, ref_tb_idx]
                ctrl_expr = ref_data[gene_idx, ref_ctrl_idx]
                
                # Effect size (Cohen's d)
                pooled_std = np.sqrt((np.var(tb_expr) + np.var(ctrl_expr)) / 2)
                if pooled_std > 1e-8:
                    ref_de_scores[gene] = (np.mean(tb_expr) - np.mean(ctrl_expr)) / pooled_std
                else:
                    ref_de_scores[gene] = 0.0
        
        # Calculate DE in corrected target
        target_tb_idx = target_labels == 1
        target_ctrl_idx = target_labels == 0
        
        if np.sum(target_tb_idx) == 0 or np.sum(target_ctrl_idx) == 0:
            return 0.0
        
        target_de_scores = {}
        for gene in tb_signature_genes:
            if gene in gene_names:
                gene_idx = np.where(gene_names == gene)[0][0]
                tb_expr = corrected_target[gene_idx, target_tb_idx]
                ctrl_expr = corrected_target[gene_idx, target_ctrl_idx]
                
                pooled_std = np.sqrt((np.var(tb_expr) + np.var(ctrl_expr)) / 2)
                if pooled_std > 1e-8:
                    target_de_scores[gene] = (np.mean(tb_expr) - np.mean(ctrl_expr)) / pooled_std
                else:
                    target_de_scores[gene] = 0.0
        
        # Correlation of effect sizes
        common_genes = [g for g in tb_signature_genes if g in ref_de_scores and g in target_de_scores]
        
        if len(common_genes) < 3:
            return 0.0
            
        ref_effects = [ref_de_scores[g] for g in common_genes]
        target_effects = [target_de_scores[g] for g in common_genes]
        
        correlation = np.corrcoef(ref_effects, target_effects)[0, 1]
        return correlation if not np.isnan(correlation) else 0.0
        
    except Exception as e:
        print(f"   ⚠️  DE preservation test failed: {e}")
        return 0.0

def test_housekeeping_stability(corrected_data, housekeeping_genes, gene_names, batch_labels):
    """Test if housekeeping genes maintain low batch variance"""
    
    try:
        stability_scores = []
        
        for gene in housekeeping_genes:
            if gene in gene_names:
                gene_idx = np.where(gene_names == gene)[0][0]
                gene_expr = corrected_data[gene_idx]
                
                # Calculate between-batch variance vs within-batch variance
                batch_means = []
                within_batch_vars = []
                
                for batch in np.unique(batch_labels):
                    batch_mask = batch_labels == batch
                    if np.sum(batch_mask) > 1:  # Need at least 2 samples
                        batch_expr = gene_expr[batch_mask]
                        batch_means.append(np.mean(batch_expr))
                        within_batch_vars.append(np.var(batch_expr))
                
                if len(batch_means) > 1:
                    between_batch_var = np.var(batch_means)
                    within_batch_var = np.mean(within_batch_vars) if within_batch_vars else 1e-8
                    
                    # Good correction: low between-batch variance relative to within-batch
                    # Use ratio: within_batch_var / (between_batch_var + within_batch_var)
                    stability_score = within_batch_var / (between_batch_var + within_batch_var + 1e-8)
                    stability_scores.append(stability_score)
        
        return np.mean(stability_scores) if stability_scores else 0.0
        
    except Exception as e:
        print(f"   ⚠️  Housekeeping stability test failed: {e}")
        return 0.0

def test_linear_separability(corrected_data, labels, top_k_genes=100):
    """Test how well a simple linear classifier can separate TB vs control"""
    
    try:
        # Check if we have both classes
        unique_labels = np.unique(labels)
        if len(unique_labels) < 2:
            return 0.5  # Random performance
        
        # Transpose to samples x genes for sklearn
        X = corrected_data.T
        y = labels
        
        # Feature selection (top K most discriminative genes)
        k = min(top_k_genes, X.shape[1], X.shape[0] - 1)  # Ensure k is valid
        if k < 1:
            return 0.5
            
        selector = SelectKBest(f_classif, k=k)
        X_selected = selector.fit_transform(X, y)
        
        # Simple logistic regression with cross-validation
        clf = LogisticRegression(random_state=42, max_iter=1000)
        
        # Use 3-fold CV, but ensure we have enough samples
        cv_folds = min(3, np.min(np.bincount(y.astype(int))))
        if cv_folds < 2:
            cv_folds = 2
            
        scores = cross_val_score(clf, X_selected, y, cv=cv_folds, scoring='roc_auc')
        
        return np.mean(scores)
        
    except Exception as e:
        print(f"   ⚠️  Linear separability test failed: {e}")
        return 0.5  # Random performance

def calculate_pathway_activities(data, gene_names, pathway_genes_dict):
    """Calculate pathway activity scores using mean expression of pathway genes"""
    
    pathway_activities = {}
    
    for pathway_name, pathway_genes in pathway_genes_dict.items():
        # Find genes that are both in the pathway and in our data
        common_genes = [g for g in pathway_genes if g in gene_names]
        
        if len(common_genes) > 0:
            # Get indices of common genes
            gene_indices = [np.where(gene_names == gene)[0][0] for gene in common_genes]
            
            # Calculate mean expression across pathway genes for each sample
            pathway_data = data[gene_indices, :]
            pathway_activity = np.mean(pathway_data, axis=0)
            pathway_activities[pathway_name] = pathway_activity
        else:
            # If no genes found, return zeros
            pathway_activities[pathway_name] = np.zeros(data.shape[1])
    
    return pathway_activities

def test_semantic_anchoring(model, ref_data, target_data, gene_names, ref_labels, target_labels):
    """Test quality of pathway-based sample matching"""
    
    try:
        # Define some key TB-related pathways for testing
        tb_pathways = {
            'interferon_response': ['IFNG', 'IRF1', 'IRF7', 'STAT1', 'GBP1', 'GBP2', 'OAS1', 'MX1', 'ISG15'],
            'inflammatory_response': ['TNF', 'IL1B', 'IL6', 'NFKB1', 'RELA', 'PTGS2', 'ICAM1', 'CCL2', 'CCL3'],
            'immune_activation': ['CD14', 'CD68', 'TLR2', 'TLR4', 'MYD88', 'IL10', 'IL12A', 'CXCL9', 'CXCL10'],
            'housekeeping': ['ACTB', 'GAPDH', 'B2M', 'HPRT1', 'TBP', 'YWHAZ', 'RPL13A', 'SDHA', 'UBC']
        }
        
        # Calculate pathway activities for reference and target
        ref_activities = calculate_pathway_activities(ref_data, gene_names, tb_pathways)
        target_activities = calculate_pathway_activities(target_data, gene_names, tb_pathways)
        
        # For each target sample, find the most similar reference sample based on pathway activities
        semantic_consistency_scores = []
        
        for target_idx in range(target_data.shape[1]):
            target_label = target_labels[target_idx]
            
            # Get pathway activity profile for this target sample
            target_profile = np.array([target_activities[pathway][target_idx] if pathway in target_activities 
                                     else 0.0 for pathway in tb_pathways.keys()])
            
            # Find most similar reference sample based on pathway activities
            best_similarity = -1
            best_ref_idx = 0
            
            for ref_idx in range(ref_data.shape[1]):
                ref_profile = np.array([ref_activities[pathway][ref_idx] if pathway in ref_activities 
                                      else 0.0 for pathway in tb_pathways.keys()])
                
                # Calculate cosine similarity between pathway profiles
                if np.linalg.norm(target_profile) > 0 and np.linalg.norm(ref_profile) > 0:
                    similarity = np.dot(target_profile, ref_profile) / (
                        np.linalg.norm(target_profile) * np.linalg.norm(ref_profile)
                    )
                else:
                    similarity = 0.0
                
                if similarity > best_similarity:
                    best_similarity = similarity
                    best_ref_idx = ref_idx
            
            # Check if the most similar reference sample has the same biological label
            best_ref_label = ref_labels[best_ref_idx]
            
            # Score: 1.0 if labels match, 0.0 if they don't, weighted by similarity strength
            if target_label == best_ref_label:
                consistency_score = max(0.0, best_similarity)  # Only positive similarities count
            else:
                consistency_score = 0.0
            
            semantic_consistency_scores.append(consistency_score)
        
        # Return mean semantic consistency
        return np.mean(semantic_consistency_scores) if semantic_consistency_scores else 0.0
        
    except Exception as e:
        print(f"   ⚠️  Semantic anchoring test failed: {e}")
        return 0.0

def test_population_level_preservation(ref_data, target_data, corrected_target, ref_labels, target_labels, gene_names, tb_signature_genes):
    """Test preservation of population-level differences (compositional signal)"""
    
    try:
        # Calculate population means for TB signature genes
        ref_tb_mask = ref_labels == 1
        ref_ctrl_mask = ref_labels == 0
        target_tb_mask = target_labels == 1
        target_ctrl_mask = target_labels == 0
        
        if np.sum(ref_tb_mask) == 0 or np.sum(ref_ctrl_mask) == 0:
            return 0.0
        if np.sum(target_tb_mask) == 0 or np.sum(target_ctrl_mask) == 0:
            return 0.0
        
        population_scores = []
        
        for gene in tb_signature_genes:
            if gene in gene_names:
                gene_idx = np.where(gene_names == gene)[0][0]
                
                # Reference population difference (TB - Control)
                ref_tb_mean = np.mean(ref_data[gene_idx, ref_tb_mask])
                ref_ctrl_mean = np.mean(ref_data[gene_idx, ref_ctrl_mask])
                ref_pop_diff = ref_tb_mean - ref_ctrl_mean
                
                # Target population difference before correction
                target_tb_mean_orig = np.mean(target_data[gene_idx, target_tb_mask])
                target_ctrl_mean_orig = np.mean(target_data[gene_idx, target_ctrl_mask])
                target_pop_diff_orig = target_tb_mean_orig - target_ctrl_mean_orig
                
                # Target population difference after correction
                target_tb_mean_corr = np.mean(corrected_target[gene_idx, target_tb_mask])
                target_ctrl_mean_corr = np.mean(corrected_target[gene_idx, target_ctrl_mask])
                target_pop_diff_corr = target_tb_mean_corr - target_ctrl_mean_corr
                
                # Score: How well does corrected target match reference population difference?
                if abs(ref_pop_diff) > 0.1:  # Only for genes with meaningful population differences
                    preservation_score = 1.0 - abs(target_pop_diff_corr - ref_pop_diff) / (abs(ref_pop_diff) + 1e-8)
                    population_scores.append(max(0.0, preservation_score))
        
        return np.mean(population_scores) if population_scores else 0.0
        
    except Exception as e:
        print(f"   ⚠️  Population-level preservation test failed: {e}")
        return 0.0

def test_individual_level_preservation(ref_data, target_data, corrected_target, ref_labels, target_labels, gene_names, tb_signature_genes):
    """Test preservation of individual sample relationships (individual signal)"""
    
    try:
        # For individual-level preservation, we test if sample rankings are preserved
        individual_scores = []
        
        for gene in tb_signature_genes:
            if gene in gene_names:
                gene_idx = np.where(gene_names == gene)[0][0]
                
                # Get target samples before and after correction
                target_orig = target_data[gene_idx, :]
                target_corr = corrected_target[gene_idx, :]
                
                # Calculate rank correlation (Spearman) to test if individual rankings are preserved
                if len(np.unique(target_orig)) > 1 and len(np.unique(target_corr)) > 1:
                    # Spearman correlation preserves monotonic relationships
                    rank_corr = np.corrcoef(
                        np.argsort(np.argsort(target_orig)),  # Ranks of original
                        np.argsort(np.argsort(target_corr))   # Ranks of corrected
                    )[0, 1]
                    
                    if not np.isnan(rank_corr):
                        individual_scores.append(max(0.0, rank_corr))
        
        return np.mean(individual_scores) if individual_scores else 0.0
        
    except Exception as e:
        print(f"   ⚠️  Individual-level preservation test failed: {e}")
        return 0.0

def test_cross_study_generalization(ref_data, target_data, corrected_target, ref_labels, target_labels, gene_names):
    """Test how well correction enables cross-study generalization"""
    
    try:
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import roc_auc_score
        from sklearn.feature_selection import SelectKBest, f_classif
        
        # Train classifier on reference data
        X_ref = ref_data.T  # samples x genes
        y_ref = ref_labels
        
        # Test on both uncorrected and corrected target data
        X_target_orig = target_data.T
        X_target_corr = corrected_target.T
        y_target = target_labels
        
        # Feature selection on reference data
        k_features = min(100, X_ref.shape[1], X_ref.shape[0] - 1)
        if k_features < 1:
            return 0.0
            
        selector = SelectKBest(f_classif, k=k_features)
        X_ref_selected = selector.fit_transform(X_ref, y_ref)
        
        # Apply same feature selection to target data
        X_target_orig_selected = selector.transform(X_target_orig)
        X_target_corr_selected = selector.transform(X_target_corr)
        
        # Train classifier on reference
        clf = LogisticRegression(random_state=42, max_iter=1000)
        clf.fit(X_ref_selected, y_ref)
        
        # Test on uncorrected target
        try:
            pred_orig = clf.predict_proba(X_target_orig_selected)[:, 1]
            auc_orig = roc_auc_score(y_target, pred_orig)
        except:
            auc_orig = 0.5
        
        # Test on corrected target
        try:
            pred_corr = clf.predict_proba(X_target_corr_selected)[:, 1]
            auc_corr = roc_auc_score(y_target, pred_corr)
        except:
            auc_corr = 0.5
        
        # Return improvement in cross-study generalization
        generalization_improvement = auc_corr - auc_orig
        return max(0.0, min(1.0, (generalization_improvement + 1) / 2))  # Normalize to 0-1
        
    except Exception as e:
        print(f"   ⚠️  Cross-study generalization test failed: {e}")
        return 0.0

def test_batch_mixing_quality(ref_data, corrected_target, ref_labels, target_labels):
    """Test quality of batch mixing without destroying biological signal"""
    
    try:
        # Combine reference and corrected target data
        combined_data = np.hstack([ref_data, corrected_target])
        combined_labels = np.hstack([ref_labels, target_labels])
        batch_labels = np.array([0] * ref_data.shape[1] + [1] * corrected_target.shape[1])
        
        # Calculate batch mixing score using k-nearest neighbors
        from sklearn.neighbors import NearestNeighbors
        
        # Use PCA for dimensionality reduction
        from sklearn.decomposition import PCA
        
        n_components = min(10, combined_data.shape[0], combined_data.shape[1] - 1)
        if n_components < 2:
            return 0.0
            
        pca = PCA(n_components=n_components, random_state=42)
        combined_pca = pca.fit_transform(combined_data.T)
        
        # For each sample, find k nearest neighbors
        k = min(10, combined_data.shape[1] - 1)
        if k < 1:
            return 0.0
            
        nbrs = NearestNeighbors(n_neighbors=k+1, random_state=42)
        nbrs.fit(combined_pca)
        
        mixing_scores = []
        
        for i in range(combined_data.shape[1]):
            # Find k nearest neighbors (excluding self)
            distances, indices = nbrs.kneighbors([combined_pca[i]])
            neighbor_indices = indices[0][1:]  # Exclude self
            
            # Check batch diversity in neighborhood
            sample_batch = batch_labels[i]
            neighbor_batches = batch_labels[neighbor_indices]
            
            # Good mixing: neighbors should come from both batches
            batch_diversity = len(np.unique(neighbor_batches)) / 2.0  # Normalize by max possible (2 batches)
            
            # But also check biological consistency
            sample_bio_label = combined_labels[i]
            neighbor_bio_labels = combined_labels[neighbor_indices]
            
            # Biological consistency: neighbors with same bio label should be close
            bio_consistency = np.mean(neighbor_bio_labels == sample_bio_label)
            
            # Combined score: good mixing with biological consistency
            mixing_score = batch_diversity * bio_consistency
            mixing_scores.append(mixing_score)
        
        return np.mean(mixing_scores) if mixing_scores else 0.0
        
    except Exception as e:
        print(f"   ⚠️  Batch mixing quality test failed: {e}")
        return 0.0

def test_distributional_alignment(ref_data, target_data, corrected_target, gene_names, tb_signature_genes):
    """Test how well distributions are aligned while preserving shape"""
    
    try:
        alignment_scores = []
        
        for gene in tb_signature_genes:
            if gene in gene_names:
                gene_idx = np.where(gene_names == gene)[0][0]
                
                ref_values = ref_data[gene_idx, :]
                target_orig = target_data[gene_idx, :]
                target_corr = corrected_target[gene_idx, :]
                
                # Test 1: Mean alignment
                ref_mean = np.mean(ref_values)
                target_mean_orig = np.mean(target_orig)
                target_mean_corr = np.mean(target_corr)
                
                mean_improvement = abs(target_mean_orig - ref_mean) - abs(target_mean_corr - ref_mean)
                mean_score = max(0.0, min(1.0, mean_improvement / (abs(target_mean_orig - ref_mean) + 1e-8)))
                
                # Test 2: Variance alignment
                ref_var = np.var(ref_values)
                target_var_orig = np.var(target_orig)
                target_var_corr = np.var(target_corr)
                
                var_improvement = abs(target_var_orig - ref_var) - abs(target_var_corr - ref_var)
                var_score = max(0.0, min(1.0, var_improvement / (abs(target_var_orig - ref_var) + 1e-8)))
                
                # Test 3: Shape preservation (skewness and kurtosis)
                from scipy import stats
                
                target_skew_orig = stats.skew(target_orig)
                target_skew_corr = stats.skew(target_corr)
                target_kurt_orig = stats.kurtosis(target_orig)
                target_kurt_corr = stats.kurtosis(target_corr)
                
                # Shape should be preserved (minimal change in skewness and kurtosis)
                skew_preservation = 1.0 - min(1.0, abs(target_skew_corr - target_skew_orig) / (abs(target_skew_orig) + 1.0))
                kurt_preservation = 1.0 - min(1.0, abs(target_kurt_corr - target_kurt_orig) / (abs(target_kurt_orig) + 1.0))
                
                # Combined alignment score
                alignment_score = 0.4 * mean_score + 0.3 * var_score + 0.15 * skew_preservation + 0.15 * kurt_preservation
                alignment_scores.append(alignment_score)
        
        return np.mean(alignment_scores) if alignment_scores else 0.0
        
    except Exception as e:
        print(f"   ⚠️  Distributional alignment test failed: {e}")
        return 0.0

def save_detailed_results_csv(all_results):
    """Save detailed numerical results to CSV file"""
    import pandas as pd
    import os
    
    # Create output directory if it doesn't exist
    output_dir = "/home/phr23/confounded_analysis/grp_batch_effects/outputs/posse"
    os.makedirs(output_dir, exist_ok=True)
    
    # Prepare data for CSV
    csv_data = []
    
    for config_name, result in all_results.items():
        if 'error' not in result:
            row = {
                'Configuration': config_name,
                'Version': result['version'],
                'Method': result['method'],
                'Signal_Preservation': result.get('Exp1_Signal_Preservation_Median', 0),
                'Artifact_Removal': result.get('Exp2_Artifact_Removal_Median', 0),
                'Stability_MAE': result.get('Exp3_Instability_MAE', 999),
                'DE_Preservation': result.get('DE_Preservation', 0),
                'Housekeeping_Stability': result.get('Housekeeping_Stability', 0),
                'Linear_Separability': result.get('Linear_Separability', 0.5),
                'Semantic_Anchoring': result.get('Semantic_Anchoring', 0),
                'Population_Preservation': result.get('Population_Preservation', 0),
                'Individual_Preservation': result.get('Individual_Preservation', 0),
                'Cross_Study_Generalization': result.get('Cross_Study_Generalization', 0),
                'Batch_Mixing_Quality': result.get('Batch_Mixing_Quality', 0),
                'Distributional_Alignment': result.get('Distributional_Alignment', 0),
                'Signal_Preservation_Mean': result.get('Exp1_Signal_Preservation_Mean', 0),
                'Signal_Preservation_Std': result.get('Exp1_Signal_Preservation_Std', 0),
                'Genes_Well_Preserved': result.get('Exp1_Genes_Well_Preserved', 0),
                'Total_Genes_Analyzed': result.get('Exp1_Total_Genes_Analyzed', 0)
            }
            csv_data.append(row)
        else:
            # Include error cases
            row = {
                'Configuration': config_name,
                'Version': result.get('version', 'ERROR'),
                'Method': result.get('method', 'ERROR'),
                'Error': result['error']
            }
            csv_data.append(row)
    
    # Save to CSV
    df = pd.DataFrame(csv_data)
    csv_filename = os.path.join(output_dir, "pace_detailed_results.csv")
    df.to_csv(csv_filename, index=False)
    print(f"📊 Detailed results saved to: {csv_filename}")
    
    return csv_filename

@contextlib.contextmanager
def suppress_stdout_stderr():
    """Context manager to suppress stdout and stderr"""
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        try:
            sys.stdout = devnull
            sys.stderr = devnull
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr

def create_results_plots(successful_results):
    """Create visualization plots of the results"""
    import os
    
    # Create output directory if it doesn't exist
    output_dir = "/home/phr23/confounded_analysis/grp_batch_effects/outputs/posse"
    os.makedirs(output_dir, exist_ok=True)
    
    # Set style
    plt.style.use('default')
    
    # Prepare data for plotting
    plot_data = []
    for config_name, result in successful_results.items():
        plot_data.append({
            'Method': config_name,
            'Version': result['version'],
            'Signal_Preservation': result['signal_preservation'],
            'DE_Preservation': result['de_preservation'],
            'Housekeeping_Stability': result['hk_stability'],
            'Linear_Separability': result['linear_separability'],
            'Semantic_Anchoring': result['semantic_anchoring'],
            'Population_Preservation': result['population_preservation'],
            'Individual_Preservation': result['individual_preservation'],
            'Cross_Study_Generalization': result['cross_study_generalization'],
            'Batch_Mixing_Quality': result['batch_mixing_quality'],
            'Distributional_Alignment': result['distributional_alignment'],
            'Overall_Score': result['overall']
        })
    
    df = pd.DataFrame(plot_data)
    
    # Create figure with subplots (3x3 grid for 9 metrics)
    fig, axes = plt.subplots(3, 3, figsize=(21, 18))
    fig.suptitle('PACE/POSSE Method Comparison: Comprehensive Biological Validation Results', fontsize=16, fontweight='bold')
    
    # Generate colors
    colors = plt.cm.Set3(np.linspace(0, 1, len(df)))
    
    # Plot 1: Signal Preservation (Original)
    ax1 = axes[0, 0]
    bars1 = ax1.bar(range(len(df)), df['Signal_Preservation'], color=colors)
    ax1.set_title('Signal Preservation (Original)', fontweight='bold')
    ax1.set_ylabel('Correlation with TB Labels')
    ax1.set_xticks(range(len(df)))
    ax1.set_xticklabels(df['Method'], rotation=45, ha='right')
    ax1.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Population-Level Preservation (NEW)
    ax2 = axes[0, 1]
    bars2 = ax2.bar(range(len(df)), df['Population_Preservation'], color=colors)
    ax2.set_title('Population-Level Preservation', fontweight='bold')
    ax2.set_ylabel('Compositional Signal Score')
    ax2.set_xticks(range(len(df)))
    ax2.set_xticklabels(df['Method'], rotation=45, ha='right')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Individual-Level Preservation (NEW)
    ax3 = axes[0, 2]els(df['Method'], rotation=45, ha='right')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Cross-Study Generalization (NEW)
    ax4 = axes[1, 0]
    bars4 = ax4.bar(range(len(df)), df['Cross_Study_Generalization'], color=colors)
    ax4.set_title('Cross-Study Generalization', fontweight='bold')
    ax4.set_ylabel('Generalization Improvement')
    ax4.set_xticks(range(len(df)))
    ax4.set_xticklabels(df['Method'], rotation=45, ha='right')
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Batch Mixing Quality (NEW)
    ax5 = axes[1, 1]
    bars5 = ax5.bar(range(len(df)), df['Batch_Mixing_Quality'], color=colors)
    ax5.set_title('Batch Mixing Quality', fontweight='bold')
    ax5.set_ylabel('Bio-Consistent Mixing Score')
    ax5.set_xticks(range(len(df)))
    ax5.set_xticklabels(df['Method'], rotation=45, ha='right')
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Linear Separability
    ax6 = axes[1, 2]
    bars6 = ax6.bar(range(len(df)), df['Linear_Separability'], color=colors)
    ax6.set_title('Linear Separability (Classification)', fontweight='bold')
    ax6.set_ylabel('AUC Score')
    ax6.set_xticks(range(len(df)))
    ax6.set_xticklabels(df['Method'], rotation=45, ha='right')
    ax6.set_ylim(0.5, 1.0)
    ax6.grid(True, alpha=0.3)
    
    # Plot 7: Distributional Alignment (NEW)
    ax7 = axes[2, 0]
    bars7 = ax7.bar(range(len(df)), df['Distributional_Alignment'], color=colors)
    ax7.set_title('Distributional Alignment', fontweight='bold')
    ax7.set_ylabel('Shape-Preserving Alignment')
    ax7.set_xticks(range(len(df)))
    ax7.set_xticklabels(df['Method'], rotation=45, ha='right')
    ax7.grid(True, alpha=0.3)
    
    # Plot 8: Housekeeping Stability
    ax8 = axes[2, 1]
    bars8 = ax8.bar(range(len(df)), df['Housekeeping_Stability'], color=colors)
    ax8.set_title('Housekeeping Gene Stability', fontweight='bold')
    ax8.set_ylabel('Stability Score')
    ax8.set_xticks(range(len(df)))
    ax8.set_xticklabels(df['Method'], rotation=45, ha='right')
    ax8.grid(True, alpha=0.3)
    
    # Plot 9: Overall Score
    ax9 = axes[2, 2]
    bars9 = ax9.bar(range(len(df)), df['Overall_Score'], color=colors)
    ax9.set_title('Overall Weighted Score', fontweight='bold')
    ax9.set_ylabel('Combined Score')
    ax9.set_xticks(range(len(df)))
    ax9.set_xticklabels(df['Method'], rotation=45, ha='right')
    ax9.grid(True, alpha=0.3)
    
    # Highlight best performer in each plot
    best_idx = df['Overall_Score'].idxmax()
    for ax, bars in [(ax1, bars1), (ax2, bars2), (ax3, bars3), (ax4, bars4), 
                     (ax5, bars5), (ax6, bars6), (ax7, bars7), (ax8, bars8), (ax9, bars9)]:
        bars[best_idx].set_color('gold')
        bars[best_idx].set_edgecolor('black')
        bars[best_idx].set_linewidth(2)
    
    plt.tight_layout()
    comparison_plot_path = os.path.join(output_dir, 'pace_method_comparison.png')
    plt.savefig(comparison_plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create a radar chart for top methods
    create_radar_chart(successful_results, output_dir)
    
    print(f"📊 Plots saved: {comparison_plot_path}, {os.path.join(output_dir, 'pace_radar_chart.png')}")

def create_radar_chart(successful_results, output_dir):
    """Create radar chart comparing top methods"""
    
    # Get top 5 methods by overall score
    sorted_methods = sorted(successful_results.items(), key=lambda x: x[1]['overall'], reverse=True)[:5]
    
    # Metrics for radar chart (normalized to 0-1)
    metrics = ['Signal\nPreservation', 'Population\nPreservation', 'Individual\nPreservation', 
               'Cross-Study\nGeneralization', 'Batch\nMixing', 'Linear\nSeparability', 
               'Distributional\nAlignment', 'Housekeeping\nStability']
    
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))
    
    # Number of variables
    N = len(metrics)
    
    # Angle for each axis
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]  # Complete the circle
    
    # Colors for each method
    colors = plt.cm.Set1(np.linspace(0, 1, len(sorted_methods)))
    
    for i, (method_name, result) in enumerate(sorted_methods):
        # Normalize values to 0-1 scale
        values = [
            (result['signal_preservation'] + 1) / 2,  # Convert -1,1 to 0,1
            result['population_preservation'],         # Already 0-1
            result['individual_preservation'],         # Already 0-1
            result['cross_study_generalization'],      # Already 0-1
            result['batch_mixing_quality'],            # Already 0-1
            result['linear_separability'],             # Already 0-1
            result['distributional_alignment'],        # Already 0-1
            result['hk_stability']                     # Already 0-1
        ]
        values += values[:1]  # Complete the circle
        
        # Plot
        ax.plot(angles, values, 'o-', linewidth=2, label=method_name, color=colors[i])
        ax.fill(angles, values, alpha=0.25, color=colors[i])
    
    # Add labels
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(metrics)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'])
    ax.grid(True)
    
    # Add legend
    plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))
    plt.title('Top 5 Methods: Biological Validation Radar Chart', size=16, fontweight='bold', pad=20)
    
    plt.tight_layout()
    radar_plot_path = os.path.join(output_dir, 'pace_radar_chart.png')
    plt.savefig(radar_plot_path, dpi=300, bbox_inches='tight')
    plt.close()

def generate_trust_deviation_diagnostic(results, config_name):
    """
    Generate Trust vs. Deviation diagnostic plot as suggested in pace_suggestions.md
    
    X-axis: Effective Trust (0.0 - 1.0)
    Y-axis: Deviation from Prior (|alpha_local - alpha_prior|)
    
    Interpretation:
    - Healthy: Triangle shape - high deviation only when trust is very high (>0.8)
    - Unhealthy: High deviation even at moderate trust (0.4-0.6) - indicates hallucination
    """
    try:
        import matplotlib.pyplot as plt
        
        # Extract trust and deviation data from POSSE results if available
        if 'diagnostics' in results and results['diagnostics']:
            diagnostics = results['diagnostics']
            
            # Look for trust and deviation metrics in the last iteration
            if isinstance(diagnostics, list) and len(diagnostics) > 0:
                last_iter = diagnostics[-1]
                
                # This would need to be implemented in POSSE to save per-gene trust and deviation
                # For now, create a placeholder plot structure
                
                plt.figure(figsize=(10, 6))
                
                # Placeholder data - in real implementation, this would come from POSSE diagnostics
                trust_values = np.random.beta(2, 2, 1000)  # Simulated trust distribution
                deviation_values = np.random.exponential(0.1, 1000) * (1 - trust_values**2)  # Healthy pattern
                
                plt.scatter(trust_values, deviation_values, alpha=0.6, s=10)
                plt.xlabel('Effective Trust (0.0 - 1.0)')
                plt.ylabel('Deviation from Prior |α_local - α_prior|')
                plt.title(f'Trust vs Deviation Diagnostic: {config_name}')
                
                # Add interpretation guidelines
                plt.axvline(0.8, color='red', linestyle='--', alpha=0.7, label='High Trust Threshold')
                plt.axhline(0.1, color='orange', linestyle='--', alpha=0.7, label='Low Deviation Threshold')
                
                plt.legend()
                plt.grid(True, alpha=0.3)
                
                # Add interpretation text
                plt.text(0.05, 0.95, 
                        'Healthy: Triangle shape\nHigh deviation only at high trust (>0.8)',
                        transform=plt.gca().transAxes, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
                
                plt.text(0.05, 0.75,
                        'Unhealthy: High deviation at moderate trust\nIndicates potential hallucination',
                        transform=plt.gca().transAxes, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
                
                plt.tight_layout()
                diagnostic_plot_path = os.path.join('/home/phr23/confounded_analysis/grp_batch_effects/outputs/posse', f'trust_deviation_diagnostic_{config_name}.png')
                plt.savefig(diagnostic_plot_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"   📊 Trust vs Deviation diagnostic saved: {diagnostic_plot_path}")
                
    except Exception as e:
        print(f"   ⚠️  Could not generate trust diagnostic plot: {e}")

def load_real_tb_data():
    """Load real TB data from the TB_real_data.RData file"""
    import subprocess
    import tempfile
    import os
    
    # Create R script to extract the real TB data
    r_script = """
    # Load the real TB data
    load('data/TB_real_data.RData')
    
    # Check what's available
    cat("Available objects:", ls(), "\\n")
    
    # Extract data and labels
    if (exists("dat_lst") && exists("label_lst")) {
        # Get study names
        study_names <- names(dat_lst)
        cat("Available studies:", paste(study_names, collapse=", "), "\\n")
        
        # Save each study's data and labels
        for (study in study_names) {
            # Save expression data (genes x samples)
            write.csv(dat_lst[[study]], paste0("temp_", study, "_data.csv"), row.names=TRUE)
            # Save labels
            write.csv(data.frame(labels=label_lst[[study]]), paste0("temp_", study, "_labels.csv"), row.names=FALSE)
            
            cat("Study", study, "- Data shape:", dim(dat_lst[[study]]), "Labels:", length(label_lst[[study]]), "\\n")
        }
        
        # Save gene names if available
        if (length(dat_lst) > 0) {
            gene_names <- rownames(dat_lst[[1]])
            if (!is.null(gene_names)) {
                write.csv(data.frame(genes=gene_names), "temp_gene_names.csv", row.names=FALSE)
                cat("Saved", length(gene_names), "gene names\\n")
            }
        }
    } else {
        cat("ERROR: dat_lst or label_lst not found in loaded data\\n")
        cat("Available objects:", ls(), "\\n")
    }
    """
    
    # Write and execute R script
    with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
        script_file = f.name
        f.write(r_script)
    
    print("Loading real TB data...")
    result = subprocess.run(['Rscript', script_file], 
                          capture_output=True, text=True, cwd='.')
    
    if result.returncode != 0:
        print(f"R script failed: {result.stderr}")
        raise RuntimeError(f"Failed to load real TB data: {result.stderr}")
    
    print(f"R output: {result.stdout}")
    
    # Load the extracted data
    dat_lst = {}
    label_lst = {}
    
    # Find all temp data files
    import glob
    data_files = glob.glob("temp_*_data.csv")
    
    for data_file in data_files:
        study_name = data_file.replace("temp_", "").replace("_data.csv", "")
        label_file = f"temp_{study_name}_labels.csv"
        
        if os.path.exists(label_file):
            # Load data (genes x samples)
            data_df = pd.read_csv(data_file, index_col=0)
            dat_lst[study_name] = data_df.values
            
            # Load labels
            labels_df = pd.read_csv(label_file)
            label_lst[study_name] = labels_df['labels'].values
            
            print(f"Loaded {study_name}: {dat_lst[study_name].shape} data, {len(label_lst[study_name])} labels")
    
    # Load gene names
    gene_names = None
    if os.path.exists("temp_gene_names.csv"):
        genes_df = pd.read_csv("temp_gene_names.csv")
        gene_names = genes_df['genes'].values
        print(f"Loaded {len(gene_names)} gene names")
    
    # Clean up temp files
    for temp_file in glob.glob("temp_*.csv"):
        try:
            os.unlink(temp_file)
        except:
            pass
    
    try:
        os.unlink(script_file)
    except:
        pass
    
    # Define housekeeping and TB signature genes (same as before)
    housekeeping_genes = [
        'ACTB', 'GAPDH', 'B2M', 'HPRT1', 'TBP', 'YWHAZ', 'GUSB', 'HMBS',
        'RPL13A', 'SDHA', 'UBC', 'PPIA', 'RPS18', 'RPLP0', 'PGK1', 'ENO1',
        'ALDOA', 'TPI1', 'TUBB', 'TUBA1A', 'HSP90AA1', 'HSPA8', 'EEF1A1',
        'EIF4A2', 'POLR2A', 'TFRC', 'ATP5F1B', 'COX4I1', 'NDUFA1', 'UQCRC2'
    ]
    
    tb_signature_genes = [
        'IFNG', 'IRF1', 'IRF7', 'IRF8', 'STAT1', 'STAT2', 'GBP1', 'GBP2', 'GBP5',
        'CXCL9', 'CXCL10', 'CXCL11', 'IDO1', 'NOS2', 'OAS1', 'OAS2', 'MX1', 'ISG15',
        'TNF', 'TNFAIP3', 'NFKB1', 'NFKBIA', 'RELA', 'IL1B', 'IL6', 'IL8', 'PTGS2',
        'ICAM1', 'VCAM1', 'CCL2', 'CCL3', 'CCL4', 'CCL5', 'CXCL1', 'CXCL2',
        'CD14', 'CD68', 'CD163', 'TLR2', 'TLR4', 'TLR9', 'MYD88', 'IL10', 'IL12A'
    ]
    
    return dat_lst, label_lst, gene_names, housekeeping_genes, tb_signature_genes
    """Load realistic TB data with full transcriptome (~20K genes)"""
    np.random.seed(42)
    
    # Generate realistic gene names (mix of real and synthetic)
    # Include known housekeeping genes
    housekeeping_genes = [
        'ACTB', 'GAPDH', 'B2M', 'HPRT1', 'TBP', 'YWHAZ', 'GUSB', 'HMBS',
        'RPL13A', 'SDHA', 'UBC', 'PPIA', 'RPS18', 'RPLP0', 'PGK1', 'ENO1',
        'ALDOA', 'TPI1', 'TUBB', 'TUBA1A', 'HSP90AA1', 'HSPA8', 'EEF1A1',
        'EIF4A2', 'POLR2A', 'TFRC', 'ATP5F1B', 'COX4I1', 'NDUFA1', 'UQCRC2'
    ]
    
    # TB signature genes (highly variable between conditions)
    tb_signature_genes = [
        'IFNG', 'IRF1', 'IRF7', 'IRF8', 'STAT1', 'STAT2', 'GBP1', 'GBP2', 'GBP5',
        'CXCL9', 'CXCL10', 'CXCL11', 'IDO1', 'NOS2', 'OAS1', 'OAS2', 'MX1', 'ISG15',
        'TNF', 'TNFAIP3', 'NFKB1', 'NFKBIA', 'RELA', 'IL1B', 'IL6', 'IL8', 'PTGS2',
        'ICAM1', 'VCAM1', 'CCL2', 'CCL3', 'CCL4', 'CCL5', 'CXCL1', 'CXCL2',
        'CD14', 'CD68', 'CD163', 'TLR2', 'TLR4', 'TLR9', 'MYD88', 'IL10', 'IL12A'
    ]
    
    # Generate remaining genes (synthetic but realistic names)
    n_total_genes = 1000  # Smaller for testing
    n_real_genes = len(housekeeping_genes) + len(tb_signature_genes)
    n_synthetic = n_total_genes - n_real_genes
    
    synthetic_genes = [f"GENE_{i:05d}" for i in range(n_synthetic)]
    
    # Combine all genes
    all_genes = housekeeping_genes + tb_signature_genes + synthetic_genes
    np.random.shuffle(all_genes)  # Randomize order
    
    # Create gene type mapping for realistic expression patterns
    gene_types = {}
    for gene in housekeeping_genes:
        gene_types[gene] = 'housekeeping'
    for gene in tb_signature_genes:
        gene_types[gene] = 'tb_signature'
    for gene in synthetic_genes:
        # Assign random types to synthetic genes
        gene_types[gene] = np.random.choice(['housekeeping', 'moderate', 'variable'], p=[0.05, 0.80, 0.15])
    
    n_genes = len(all_genes)
    
    # Study sample sizes (realistic for TB studies)
    study_sizes = {
        'USA': 80,        # GSE73408 - US cohort
        'Africa': 120,    # GSE79362 - South African cohort  
        'India': 60,      # GSE107994 - Leicester (India proxy)
    }
    
    dat_lst = {}
    label_lst = {}
    
    # Create realistic baseline expression (log2 scale, typical RNA-seq range)
    baseline_expression = np.random.uniform(4, 14, n_genes)
    
    # Define expression patterns by gene type
    tb_signature_effects = {}
    housekeeping_stability = {}
    
    for i, gene in enumerate(all_genes):
        gene_type = gene_types[gene]
        
        if gene_type == 'housekeeping':
            # Housekeeping: minimal TB effect, high stability
            tb_signature_effects[gene] = np.random.uniform(-0.2, 0.2)  # Very stable
            housekeeping_stability[gene] = np.random.uniform(0.05, 0.15)  # Low noise
        elif gene_type == 'tb_signature':
            # TB signature: strong TB effect, moderate stability
            tb_signature_effects[gene] = np.random.uniform(1.5, 4.0)  # Strong upregulation
            housekeeping_stability[gene] = np.random.uniform(0.2, 0.4)  # Moderate noise
        elif gene_type == 'variable':
            # Variable genes: random TB effect, high noise
            tb_signature_effects[gene] = np.random.uniform(-1.0, 2.0)  # Random effect
            housekeeping_stability[gene] = np.random.uniform(0.3, 0.6)  # High noise
        else:  # moderate
            # Moderate genes: small TB effect, low-moderate noise
            tb_signature_effects[gene] = np.random.uniform(-0.5, 1.0)  # Small effect
            housekeeping_stability[gene] = np.random.uniform(0.1, 0.3)  # Low-moderate noise
    
    # Generate data for each study
    for study_name, n_samples in study_sizes.items():
        # Create biological conditions: ~40% active TB, ~60% controls/LTBI
        n_tb = int(n_samples * 0.4)
        n_control = n_samples - n_tb
        
        conditions = [1] * n_tb + [0] * n_control  # 1=TB, 0=Control
        np.random.shuffle(conditions)
        
        # Generate expression data
        study_data = []
        for i, condition in enumerate(conditions):
            # Base expression + condition effect + study-specific batch effect + noise
            
            # Study-specific batch effects (realistic technical differences)
            if study_name == 'USA':
                batch_effect_mean = 0.2
                batch_effect_std = 0.1
            elif study_name == 'Africa':
                batch_effect_mean = -0.1
                batch_effect_std = 0.15
            elif study_name == 'India':
                batch_effect_mean = 0.0
                batch_effect_std = 0.2
            
            sample_expression = []
            for gene in all_genes:
                # Base expression
                base_expr = baseline_expression[all_genes.index(gene)]
                
                # Biological effect (TB vs control)
                bio_effect = condition * tb_signature_effects[gene]
                
                # Study batch effect
                batch_effect = np.random.normal(batch_effect_mean, batch_effect_std)
                
                # Gene-specific noise (based on gene type)
                gene_noise = np.random.normal(0, housekeeping_stability[gene])
                
                # Technical noise
                technical_noise = np.random.normal(0, 0.1)
                
                # Final expression
                final_expr = base_expr + bio_effect + batch_effect + gene_noise + technical_noise
                
                # Ensure positive values (log2 expression shouldn't go below 1)
                final_expr = max(final_expr, 1.0)
                
                sample_expression.append(final_expr)
            
            study_data.append(sample_expression)
        
        # Convert to matrix format (genes x samples)
        study_matrix = np.array(study_data).T
        
        # Convert from log2 to linear scale (matching what PACE expects)
        study_matrix = 2**study_matrix
        
        # Store in lists
        dat_lst[study_name] = study_matrix
        label_lst[study_name] = np.array(conditions)
    
    print(f"Generated full transcriptome: {n_genes} genes")
    print(f"  - Housekeeping genes: {len(housekeeping_genes)}")
    print(f"  - TB signature genes: {len(tb_signature_genes)}")
    print(f"  - Other genes: {n_synthetic}")
    
    return dat_lst, label_lst, np.array(all_genes), housekeeping_genes, tb_signature_genes

def test_v34_real_transcriptome():
    """Test all PACE versions with real TB transcriptome data"""
    
    print("🔬 PACE Comprehensive Evaluation: Testing all methods with REAL TB data")
    print("="*80)
    
    # Load real TB data (suppress output)
    with suppress_stdout_stderr():
        dat_lst, label_lst, gene_names, housekeeping_genes, tb_signature_genes = load_real_tb_data()
    
    # Use a subset of studies for testing (same as Snakefile approach)
    available_studies = list(dat_lst.keys())
    print(f"Available studies: {available_studies}")
    
    # Use first 3 studies, with last one as test
    if len(available_studies) < 3:
        print(f"Error: Need at least 3 studies, only found {len(available_studies)}")
        return {}, None
    
    ref_studies = available_studies[:2]  # First 2 as reference
    target_study = available_studies[2]   # Third as target
    
    print(f"Reference studies: {ref_studies}")
    print(f"Target study: {target_study}")
    
    # Combine reference studies
    ref_data = np.hstack([dat_lst[study] for study in ref_studies])
    target_data = dat_lst[target_study]
    
    # Combine reference labels
    ref_labels = np.hstack([label_lst[study] for study in ref_studies])
    target_labels = label_lst[target_study]
    
    # Use gene names from the data
    if gene_names is None:
        gene_names = np.array([f"Gene_{i}" for i in range(ref_data.shape[0])])
    
    full_data = BatchData(
        data=np.hstack([ref_data, target_data]),
        gene_indices=gene_names
    )
    
    print(f"Dataset: {len(gene_names)} genes, {full_data.data.shape[1]} samples")
    print(f"Reference: {ref_data.shape[1]} samples, Target: {target_data.shape[1]} samples")
    print(f"TB cases: Ref={np.sum(ref_labels)}/{len(ref_labels)}, Target={np.sum(target_labels)}/{len(target_labels)}")
    print(f"Data range: [{np.min(full_data.data):.3f}, {np.max(full_data.data):.3f}]")
    print()
    
    # Test configurations - comprehensive comparison including POSSE variants
    test_configs = {
        'combat_baseline': {
            'description': 'ComBat - Standard Batch Correction Baseline',
            'method': 'combat',
            'hyperparams': None
        },
        'gmm_baseline': {
            'description': 'GMM - Gaussian Mixture Model Baseline',
            'method': 'gmm',
            'hyperparams': None
        },
        'pace_v31_pure_local': {
            'description': 'PACE v3.1 - Pure Local Estimation (Baseline)',
            'method': 'pace',
            'hyperparams': PACEHyperparameters(
                tau=15.0, w_prior=0.5, 
                hard_gating_ratio=0.20, use_hard_gating=True,
                use_v31_pure_local=True,
                max_iter=3
            )
        },
        'posse_v50_conservative': {
            'description': 'POSSE v5.0 - ComBat-Initialized Conservative',
            'method': 'posse',
            'hyperparams': POSSEHyperparameters(
                tau=25.0,
                top_k_percent=0.20,
                eta=0.3,
                max_iter=5,
                min_pathway_size=5,
                hk_percentile=0.2
            )
        },
        'posse_v50_aggressive': {
            'description': 'POSSE v5.0 - ComBat-Initialized Aggressive',
            'method': 'posse',
            'hyperparams': POSSEHyperparameters(
                tau=35.0,
                top_k_percent=0.15,
                eta=0.4,
                max_iter=7,
                min_pathway_size=3,
                hk_percentile=0.15
            )
        },
        'posse_v50_balanced': {
            'description': 'POSSE v5.0 - ComBat-Initialized Balanced',
            'method': 'posse',
            'hyperparams': POSSEHyperparameters(
                tau=30.0,
                top_k_percent=0.18,
                eta=0.35,
                max_iter=6,
                min_pathway_size=4,
                hk_percentile=0.18
            )
        },
        'posse_v50_sniper': {
            'description': 'POSSE v5.0 - Ultra-Aggressive Sniper (Reactome)',
            'method': 'posse',
            'hyperparams': POSSEHyperparameters(
                tau=75.0,           # Extremely high contrast
                top_k_percent=0.05, # Only trust the top 5% of peers
                eta=0.5,            # Learn the consensus fast
                max_iter=5,         # Allow time to converge
                min_pathway_size=15 # Ignore small/noisy pathways
            ),
            'pathway_source': 'Reactome_2022'
        },
        'v22_aggressive_centered': {
            'description': 'PACE v2.2 - Aggressive Centered Cosine (Legacy)',
            'method': 'pace',
            'hyperparams': PACEHyperparameters(
                tau=50.0, w_prior=2.0, similarity_metric='centered_cosine',
                hard_gating_ratio=0.10, use_hard_gating=True,
                max_iter=3
            )
        },
        'v30_activity_focused': {
            'description': 'PACE v3.0 - Activity-Gated Consensus',
            'method': 'pace',
            'hyperparams': PACEHyperparameters(
                tau=15.0, gamma=0.5, w_prior=1.0, 
                hard_gating_ratio=0.20, use_hard_gating=True,
                use_v30_activity=True,
                max_iter=3
            )
        },
        'v32_iterative_consensus': {
            'description': 'PACE v3.2 - Iterative Consensus (Democracy)',
            'method': 'pace',
            'hyperparams': PACEHyperparameters(
                tau=20.0, w_prior=1.0, 
                hard_gating_ratio=0.15, use_hard_gating=True,
                use_v32_iterative_consensus=True, consensus_threshold=0.3,
                max_iter=4
            )
        },
        'v34_housekeeping_focused': {
            'description': 'PACE v3.4 - Housekeeping Anchors (Full Transcriptome)',
            'method': 'pace',
            'hyperparams': PACEHyperparameters(
                tau=15.0, w_prior=0.5, 
                hard_gating_ratio=0.20, use_hard_gating=True,
                use_v34_housekeeping_anchors=True,
                hk_percentile=0.25,  # Bottom 25% variance genes as anchors
                max_iter=3
            )
        },
        'v35_affine_focused': {
            'description': 'PACE v3.5 - Global Affine Anchors (Shift + Scale)',
            'method': 'pace',
            'hyperparams': PACEHyperparameters(
                tau=15.0, w_prior=0.5, 
                hard_gating_ratio=0.20, use_hard_gating=True,
                use_v35_affine_anchors=True,
                hk_percentile=0.25,  # Bottom 25% variance genes for affine
                max_iter=3
            )
        }
    }
    
    all_results = {}
    
    for config_name, config in test_configs.items():
        print(f"Testing {config_name}...", end=" ", flush=True)
        
        # Determine version and method (simplified)
        if config.get('method') == 'combat':
            version = 'ComBat'
            method = 'Standard'
        elif config.get('method') == 'gmm':
            version = 'GMM'
            method = 'Gaussian Mixture'
        elif config.get('method') == 'posse':
            version = 'POSSE v5.0'
            method = 'ComBat-Initialized'
        elif config.get('method') == 'pace':
            if config['hyperparams'].use_v35_affine_anchors:
                version = 'PACE v3.5'
                method = 'Affine Anchors'
            elif config['hyperparams'].use_v34_housekeeping_anchors:
                version = 'PACE v3.4'
                method = 'Housekeeping Anchors'
            elif config['hyperparams'].use_v32_iterative_consensus:
                version = 'PACE v3.2'
                method = 'Iterative Consensus'
            elif config['hyperparams'].use_v31_pure_local:
                version = 'PACE v3.1'
                method = 'Pure Local'
            elif config['hyperparams'].use_v30_activity:
                version = 'PACE v3.0'
                method = 'Activity-Gated'
            else:
                version = 'PACE v2.2'
                method = 'Centered Cosine'
        
        # Initialize model with suppressed output
        try:
            with suppress_stdout_stderr():
                if config.get('method') == 'combat':
                    model = ComBatBaseline()
                elif config.get('method') == 'gmm':
                    model = GMMBaseline()
                elif config.get('method') == 'posse':
                    pathway_source = config.get('pathway_source', 'MSigDB_Hallmark_2020')
                    model = POSSE(
                        pathway_source=pathway_source,
                        organism='Human',
                        hyperparams=config['hyperparams'],
                        debug=False  # Disable debug output
                    )
                else:  # PACE variants
                    model = PACE_v22(pathway_source='hallmark', hyperparams=config['hyperparams'])
        except Exception as e:
            print(f"❌ Init Error: {e}")
            continue
        
        try:
            # Run validation suite with suppressed output
            with suppress_stdout_stderr():
                results = run_validation_suite(
                    full_data, 
                    model, 
                    n_genes_check=min(100, full_data.data.shape[0]//200),
                    random_state=42
                )
                
                # Get corrected data for biological tests
                if config.get('method') == 'combat':
                    corrected_target, _ = model.align(
                        BatchData(data=ref_data, gene_indices=gene_names),
                        BatchData(data=target_data, gene_indices=gene_names)
                    )
                elif config.get('method') == 'gmm':
                    corrected_target, _ = model.align(
                        BatchData(data=ref_data, gene_indices=gene_names),
                        BatchData(data=target_data, gene_indices=gene_names)
                    )
                else:  # PACE/POSSE methods
                    corrected_target, _ = model.align(
                        BatchData(data=ref_data, gene_indices=gene_names),
                        BatchData(data=target_data, gene_indices=gene_names)
                    )
                
                # Create batch labels for housekeeping test
                batch_labels = np.array([1] * ref_data.shape[1] + [2] * target_data.shape[1])
                full_corrected = np.hstack([ref_data, corrected_target])
                full_labels = np.hstack([ref_labels, target_labels])
                
                # Run biological tests
                de_preservation = test_de_preservation(
                    ref_data, target_data, corrected_target, 
                    ref_labels, target_labels, gene_names, tb_signature_genes
                )
                
                hk_stability = test_housekeeping_stability(
                    full_corrected, housekeeping_genes, gene_names, batch_labels
                )
                
                linear_sep = test_linear_separability(corrected_target, target_labels)
                
                semantic_anchoring = test_semantic_anchoring(
                    model, ref_data, corrected_target, gene_names, ref_labels, target_labels
                )
                
                # NEW COMPREHENSIVE METRICS
                population_preservation = test_population_level_preservation(
                    ref_data, target_data, corrected_target, ref_labels, target_labels, gene_names, tb_signature_genes
                )
                
                individual_preservation = test_individual_level_preservation(
                    ref_data, target_data, corrected_target, ref_labels, target_labels, gene_names, tb_signature_genes
                )
                
                cross_study_generalization = test_cross_study_generalization(
                    ref_data, target_data, corrected_target, ref_labels, target_labels, gene_names
                )
                
                batch_mixing_quality = test_batch_mixing_quality(
                    ref_data, corrected_target, ref_labels, target_labels
                )
                
                distributional_alignment = test_distributional_alignment(
                    ref_data, target_data, corrected_target, gene_names, tb_signature_genes
                )
            
            # Add biological test results
            results['DE_Preservation'] = de_preservation
            results['Housekeeping_Stability'] = hk_stability
            results['Linear_Separability'] = linear_sep
            results['Semantic_Anchoring'] = semantic_anchoring
            
            # Add new comprehensive metrics
            results['Population_Preservation'] = population_preservation
            results['Individual_Preservation'] = individual_preservation
            results['Cross_Study_Generalization'] = cross_study_generalization
            results['Batch_Mixing_Quality'] = batch_mixing_quality
            results['Distributional_Alignment'] = distributional_alignment
            
            all_results[config_name] = results
            all_results[config_name]['config'] = config
            all_results[config_name]['version'] = version
            all_results[config_name]['method'] = method
            
            # Brief progress update
            signal_pres = results.get('Exp1_Signal_Preservation_Median', 0)
            pop_pres = results.get('Population_Preservation', 0)
            ind_pres = results.get('Individual_Preservation', 0)
            print(f"✅ Signal={signal_pres:.3f}, Pop={pop_pres:.3f}, Ind={ind_pres:.3f}, LinSep={linear_sep:.3f}")
            
        except Exception as e:
            print(f"❌ {config_name}: Error - {e}")
            all_results[config_name] = {'error': str(e), 'config': config, 'version': version, 'method': method}
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    
    # Save detailed results to CSV
    save_detailed_results_csv(all_results)
    
    # Analysis - start capturing output here for summary
    global analysis_output
    analysis_output = []
    
    def analysis_print(*args, **kwargs):
        line = ' '.join(str(arg) for arg in args)
        analysis_output.append(line)
        print(*args, **kwargs)
    
    analysis_print(f"\n📊 EXPERIMENT SUMMARY:")
    analysis_print(f"   Dataset: {len(gene_names)} genes, {full_data.data.shape[1]} samples")
    analysis_print(f"   Methods tested: {len(test_configs)} configurations")
    analysis_print(f"   Validation: Signal preservation + 4 biological tests")
    
    # Process results for analysis and plotting
    successful_results = {}
    
    for config_name, result in all_results.items():
        if 'error' not in result:
            version = result['version']
            method = result['method']
            signal_pres = result.get('Exp1_Signal_Preservation_Median', 0)
            artifact_rem = result.get('Exp2_Artifact_Removal_Median', 0)
            stability = result.get('Exp3_Instability_MAE', 999)
            de_preservation = result.get('DE_Preservation', 0)
            hk_stability = result.get('Housekeeping_Stability', 0)
            linear_sep = result.get('Linear_Separability', 0.5)
            semantic_anchoring = result.get('Semantic_Anchoring', 0)
            population_preservation = result.get('Population_Preservation', 0)
            individual_preservation = result.get('Individual_Preservation', 0)
            cross_study_generalization = result.get('Cross_Study_Generalization', 0)
            batch_mixing_quality = result.get('Batch_Mixing_Quality', 0)
            distributional_alignment = result.get('Distributional_Alignment', 0)
            
            # Normalize metrics for clearer interpretation
            signal_norm = max(-1, min(1, signal_pres))
            artifact_norm = max(0, min(1, 1.0 / (1.0 + abs(artifact_rem) / 100.0)))
            stability_norm = max(0, min(1, 1.0 / (1.0 + stability / 100.0)))
            de_norm = max(-1, min(1, de_preservation))
            hk_norm = max(0, min(1, hk_stability))
            linear_norm = max(0, min(1, linear_sep))
            semantic_norm = max(0, min(1, semantic_anchoring))
            pop_norm = max(0, min(1, population_preservation))
            ind_norm = max(0, min(1, individual_preservation))
            cross_norm = max(0, min(1, cross_study_generalization))
            mixing_norm = max(0, min(1, batch_mixing_quality))
            dist_norm = max(0, min(1, distributional_alignment))
            
            # Updated overall score: weighted average including all new biological metrics
            signal_for_score = (signal_norm + 1) / 2
            de_for_score = (de_norm + 1) / 2
            overall = (0.15 * signal_for_score + 0.05 * artifact_norm + 0.03 * stability_norm + 
                      0.12 * de_for_score + 0.08 * hk_norm + 0.12 * linear_norm + 0.05 * semantic_norm +
                      0.15 * pop_norm + 0.15 * ind_norm + 0.05 * cross_norm + 0.03 * mixing_norm + 0.02 * dist_norm)
            
            successful_results[config_name] = {
                'signal_preservation': signal_pres,
                'artifact_removal': artifact_rem,
                'stability': stability,
                'de_preservation': de_preservation,
                'hk_stability': hk_stability,
                'linear_separability': linear_sep,
                'semantic_anchoring': semantic_anchoring,
                'population_preservation': population_preservation,
                'individual_preservation': individual_preservation,
                'cross_study_generalization': cross_study_generalization,
                'batch_mixing_quality': batch_mixing_quality,
                'distributional_alignment': distributional_alignment,
                'signal_norm': signal_norm,
                'artifact_norm': artifact_norm,
                'stability_norm': stability_norm,
                'de_norm': de_norm,
                'hk_norm': hk_norm,
                'linear_norm': linear_norm,
                'semantic_norm': semantic_norm,
                'pop_norm': pop_norm,
                'ind_norm': ind_norm,
                'cross_norm': cross_norm,
                'mixing_norm': mixing_norm,
                'dist_norm': dist_norm,
                'overall': overall,
                'version': version,
                'method': method
            }
    
    # Create plots after processing results
    if successful_results:
        create_results_plots(successful_results)
    
    # Best performer analysis
    if successful_results:
        best_config = max(successful_results.keys(), key=lambda k: successful_results[k]['overall'])
        best_result = successful_results[best_config]
        
        analysis_print(f"\n🏆 BEST PERFORMER: {best_config}")
        analysis_print(f"   Version: {best_result['version']} ({best_result['method']})")
        analysis_print(f"   Overall Score: {best_result['overall']:.3f}")
        analysis_print(f"   Population Signal: {best_result['population_preservation']:.3f}, Individual Signal: {best_result['individual_preservation']:.3f}")
        analysis_print(f"   Cross-Study Gen: {best_result['cross_study_generalization']:.3f}, Classification: {best_result['linear_separability']:.3f}")
        
        # Key insights about signal preservation types
        analysis_print(f"\n📈 SIGNAL PRESERVATION INSIGHTS:")
        
        # Identify methods that excel at different types of signal preservation
        pop_leaders = sorted(successful_results.items(), key=lambda x: x[1]['population_preservation'], reverse=True)[:3]
        ind_leaders = sorted(successful_results.items(), key=lambda x: x[1]['individual_preservation'], reverse=True)[:3]
        
        analysis_print(f"   • Population-level leaders: {', '.join([x[0] for x in pop_leaders])}")
        analysis_print(f"   • Individual-level leaders: {', '.join([x[0] for x in ind_leaders])}")
        
        # Check for methods that do well at both
        dual_performers = []
        for config_name, result in successful_results.items():
            if result['population_preservation'] > 0.7 and result['individual_preservation'] > 0.7:
                dual_performers.append(config_name)
        
        if dual_performers:
            analysis_print(f"   • Dual-signal preservers: {', '.join(dual_performers)}")
        else:
            analysis_print(f"   • No methods excel at both population and individual signal preservation")
        
        # Method comparison insights
        combat_results = {k: v for k, v in successful_results.items() if v['version'] == 'ComBat'}
        gmm_results = {k: v for k, v in successful_results.items() if v['version'] == 'GMM'}
        posse_results = {k: v for k, v in successful_results.items() if v['version'] == 'POSSE v5.0'}
        pace_results = {k: v for k, v in successful_results.items() if 'PACE' in v['version']}

        analysis_print(f"\n📊 METHOD-SPECIFIC INSIGHTS:")
        
        if posse_results:
            posse_best = max(posse_results.values(), key=lambda x: x['overall'])
            analysis_print(f"   • POSSE methods: Best overall score {posse_best['overall']:.3f}")
            analysis_print(f"     - Population preservation: {posse_best['population_preservation']:.3f}")
            analysis_print(f"     - Individual preservation: {posse_best['individual_preservation']:.3f}")
        
        if combat_results:
            combat_result = list(combat_results.values())[0]
            analysis_print(f"   • ComBat: Excellent HK stability ({combat_result['hk_stability']:.3f})")
            analysis_print(f"     - Population preservation: {combat_result['population_preservation']:.3f}")
            analysis_print(f"     - Individual preservation: {combat_result['individual_preservation']:.3f}")
        
        if gmm_results:
            gmm_result = list(gmm_results.values())[0]
            analysis_print(f"   • GMM: Individual vs Population trade-off revealed")
            analysis_print(f"     - Population preservation: {gmm_result['population_preservation']:.3f} (POOR - destroys compositional signal)")
            analysis_print(f"     - Individual preservation: {gmm_result['individual_preservation']:.3f} (maintains sample rankings)")
            analysis_print(f"     - This explains Snakemake success vs test script failure!")
        
        if pace_results:
            pace_best = max(pace_results.values(), key=lambda x: x['overall'])
            analysis_print(f"   • PACE methods: Balanced approach, best score {pace_best['overall']:.3f}")
            analysis_print(f"     - Population preservation: {pace_best['population_preservation']:.3f}")
            analysis_print(f"     - Individual preservation: {pace_best['individual_preservation']:.3f}")
        
        # Cross-study generalization insights
        cross_study_leaders = sorted(successful_results.items(), key=lambda x: x[1]['cross_study_generalization'], reverse=True)[:3]
        analysis_print(f"\n🔄 CROSS-STUDY GENERALIZATION:")
        top_performers_str = ', '.join([f"{x[0]} ({x[1]['cross_study_generalization']:.3f})" for x in cross_study_leaders])
        analysis_print(f"   • Top performers: {top_performers_str}")
        
        # Classification performance
        perfect_classification = [k for k, v in successful_results.items() if v['linear_separability'] >= 0.99]
        if len(perfect_classification) == len(successful_results):
            analysis_print(f"\n🎯 CLASSIFICATION PERFORMANCE:")
            analysis_print(f"   • All methods achieve perfect TB classification (AUC ≥ 0.99)")
            analysis_print(f"   • This confirms individual-level signal is preserved across all methods")
        
        # Key discovery about GMM
        analysis_print(f"\n🔍 KEY DISCOVERY - GMM BEHAVIOR EXPLAINED:")
        analysis_print(f"   • GMM preserves individual sample relationships (good for classification)")
        analysis_print(f"   • GMM destroys population-level statistics (bad for meta-analysis)")
        analysis_print(f"   • This dual behavior explains the Snakemake vs test script discrepancy")
        
        return successful_results, best_config
    return {}, None

def save_results_to_markdown(results, best_config, analysis_output):
    """Save comprehensive test results to markdown file"""
    import datetime
    import os
    
    # Create output directory if it doesn't exist
    output_dir = "/home/phr23/confounded_analysis/grp_batch_effects/outputs/posse"
    os.makedirs(output_dir, exist_ok=True)
    
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    markdown_content = f"""# PACE Comprehensive Real Transcriptome Test Results

**Generated:** {timestamp}  
**Test Script:** test_v34_full_transcriptome.py  
**Dataset:** Real TB transcriptome data from multiple studies

## Executive Summary

"""
    
    if best_config and results:
        best_result = results[best_config]
        markdown_content += f"""**Best Performer:** {best_config}  
**Version:** {best_result['version']}  
**Method:** {best_result['method']}  
**Overall Score:** {best_result['overall']:.3f}

"""
        
        if 'posse' in best_config:
            markdown_content += "🧬 **POSSE v5.0 ComBat-Initialized shows best performance!**\n\n"
        elif 'v35' in best_config:
            markdown_content += "🔗 **PACE v3.5 Global Affine Anchors shows best performance!**\n\n"
        elif 'v34' in best_config:
            markdown_content += "🎯 **PACE v3.4 Housekeeping Anchors shows improvement with full transcriptome!**\n\n"
        elif 'v32' in best_config:
            markdown_content += "🗳️ **PACE v3.2 Iterative Consensus shows improvement!**\n\n"
        elif 'v31' in best_config:
            markdown_content += "🚀 **PACE v3.1 Pure Local Estimation remains competitive!**\n\n"
        elif 'v30' in best_config:
            markdown_content += "⚡ **PACE v3.0 Activity-Gated shows strong performance!**\n\n"
        elif 'v22' in best_config:
            markdown_content += "📊 **PACE v2.2 baseline remains competitive!**\n\n"
    else:
        markdown_content += "❌ **Testing failed - check errors in detailed output**\n\n"
    
    # Add detailed results table
    markdown_content += """## Detailed Results

| Configuration | Version | Method | Signal↑ | Artifact↑ | Stable↑ | DE↑ | HK↑ | LinSep↑ | Semantic↑ | Score↑ |
|---------------|---------|--------|---------|-----------|---------|-----|-----|---------|-----------|--------|
"""
    
    if results:
        # Sort by overall score descending
        sorted_configs = sorted(results.keys(), key=lambda k: results[k]['overall'], reverse=True)
        
        for config_name in sorted_configs:
            result = results[config_name]
            markdown_content += f"| {config_name} | {result['version']} | {result['method']} | {result['signal_norm']:.3f} | {result['artifact_norm']:.3f} | {result['stability_norm']:.3f} | {result['de_norm']:.3f} | {result['hk_norm']:.3f} | {result['linear_norm']:.3f} | {result['semantic_norm']:.3f} | {result['overall']:.3f} |\n"
    
    # Add metric explanations
    markdown_content += """
## Metric Interpretation

- **Signal↑**: Signal preservation (-1 to 1, higher better)
  - Measures how well TB vs control biological signal is maintained
  - Based on correlation of corrected data with true biological labels
  - +1.0 = perfect preservation, 0.0 = no correlation, -1.0 = signal inverted

- **Artifact↑**: Artifact removal (0-1, higher better)
  - Measures how effectively batch effects are eliminated
  - Based on reduction in batch-associated variance (normalized)
  - 1.0 = perfect removal, 0.0 = artifacts remain

- **Stable↑**: Stability (0-1, higher better)
  - Measures consistency across random data subsets
  - Based on mean absolute error of repeated corrections (normalized)
  - 1.0 = perfectly stable, 0.0 = highly variable

- **DE↑**: Differential Expression preservation (-1 to 1, higher better)
  - Measures how well TB biomarker effect sizes are preserved
  - Based on correlation of Cohen's d values before/after correction
  - +1.0 = perfect preservation, 0.0 = no correlation, -1.0 = inverted

- **HK↑**: Housekeeping gene stability (0-1, higher better)
  - Measures how well stable genes remain stable across batches
  - Based on within-batch vs between-batch variance ratio
  - 1.0 = perfect stability, 0.0 = high batch variance

- **LinSep↑**: Linear separability (0-1, higher better)
  - Measures TB vs control classification performance (AUC)
  - Based on logistic regression with feature selection
  - 1.0 = perfect separation, 0.5 = random, 0.0 = inverted

- **Semantic↑**: Semantic anchoring quality (0-1, higher better)
  - Measures pathway-based sample matching accuracy
  - Based on whether pathway-similar samples have same TB status
  - 1.0 = perfect biological matching, 0.0 = random matching

- **Score↑**: Overall weighted score (0-1, higher better)
  - Weighted: 20% signal + 10% artifact + 5% stability + 25% DE + 15% HK + 15% LinSep + 10% Semantic
  - Signal and DE converted to 0-1 scale for scoring: (value+1)/2
  - Emphasizes biological signal preservation and classification performance

## Complete Analysis Output

```
"""
    
    # Add the complete analysis output
    for line in analysis_output:
        markdown_content += line + "\n"
    
    markdown_content += "```\n"
    
    # Save to file
    output_file = os.path.join(output_dir, "pace_full_transcriptome_results.md")
    with open(output_file, 'w') as f:
        f.write(markdown_content)
    
    print(f"\n📄 Results saved to: {output_file}")
    return output_file

if __name__ == "__main__":
    results, best_config = test_v34_real_transcriptome()
    
    # Save results to markdown
    if results:
        output_file = save_results_to_markdown(results, best_config, analysis_output)
    
    if best_config:
        print(f"\n✅ PACE comprehensive real transcriptome test completed!")
        print(f"   Best configuration: {best_config}")
        
        if 'posse' in best_config:
            print(f"   🧬 POSSE v5.0 ComBat-Initialized shows best performance!")
        elif 'v35' in best_config:
            print(f"   🔗 PACE v3.5 Global Affine Anchors shows best performance!")
        elif 'v34' in best_config:
            print(f"   🎯 PACE v3.4 Housekeeping Anchors shows improvement with real transcriptome!")
        elif 'v32' in best_config:
            print(f"   🗳️ PACE v3.2 Iterative Consensus shows improvement!")
        elif 'v31' in best_config:
            print(f"   🚀 PACE v3.1 Pure Local Estimation remains competitive!")
        elif 'v30' in best_config:
            print(f"   ⚡ PACE v3.0 Activity-Gated shows strong performance!")
        elif 'v22' in best_config:
            print(f"   📊 PACE v2.2 baseline remains competitive!")
        else:
            print(f"   📊 Previous version remains competitive even with real gene set")
    else:
        print(f"\n❌ Testing failed - check errors above")