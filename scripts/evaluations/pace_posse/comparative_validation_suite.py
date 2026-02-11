#!/usr/bin/env python3
"""
Comparative Validation Suite
Tests POSSE vs ComBat vs GMM to establish realistic performance baselines
"""

import sys
sys.path.append('scripts')

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import contextlib
import os
import subprocess

# Import existing modules
from pace import BatchData
from posse import POSSE, POSSEHyperparameters
from posse_v6 import POSSEv6, POSSEv6Hyperparameters, BatchData as BatchDataV6

@dataclass
class ComparativeConfig:
    """Configuration for comparative validation"""
    output_dir: str = "comparative_validation"
    save_plots: bool = True
    verbose: bool = True
    random_state: int = 42
    # Realistic thresholds (not perfection)
    silence_rmse_threshold: float = 1e-3  # More realistic than 1e-9
    alpha_deviation_threshold: float = 0.5  # More lenient
    recovery_error_threshold: float = 0.3   # Allow some error
    topology_preservation_threshold: float = 0.1  # Bimodality preservation threshold

class SuppressOutput:
    """Context manager to suppress stdout and stderr"""
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

class ComparativeValidationSuite:
    """Comparative validation suite for batch correction methods"""
    
    def __init__(self, config: ComparativeConfig = None):
        self.config = config or ComparativeConfig()
        self.results = {}
        self.method_failures = {'POSSE': [], 'POSSEv6': [], 'ComBat': [], 'GMM': []}
        
    def run_combat_correction(self, ref_data, target_data):
        """Run ComBat correction using R"""
        try:
            # Save data for R (genes x samples format)
            np.savetxt('temp_ref_data.txt', ref_data.data, delimiter='\t')  # genes x samples
            np.savetxt('temp_target_data.txt', target_data.data, delimiter='\t')  # genes x samples
            
            # R script for ComBat
            r_script = """
            library(sva)
            
            # Load data (genes x samples)
            ref_data <- read.table('temp_ref_data.txt', sep='\\t')
            target_data <- read.table('temp_target_data.txt', sep='\\t')
            
            cat("Ref data shape:", dim(ref_data), "\\n")
            cat("Target data shape:", dim(target_data), "\\n")
            
            # Combine data (genes x samples)
            combined_data <- cbind(ref_data, target_data)
            batch <- c(rep(1, ncol(ref_data)), rep(2, ncol(target_data)))
            
            cat("Combined data shape:", dim(combined_data), "\\n")
            cat("Batch vector length:", length(batch), "\\n")
            
            # Run ComBat with output suppression
            sink('/dev/null')
            corrected_data <- ComBat(dat=as.matrix(combined_data), batch=batch)
            sink()
            
            # Extract corrected target data (genes x samples)
            n_ref_samples <- ncol(ref_data)
            corrected_target <- corrected_data[, (n_ref_samples+1):ncol(corrected_data)]
            
            # Save results
            write.table(corrected_target, 'temp_combat_result.txt', sep='\\t', 
                       row.names=FALSE, col.names=FALSE)
            
            # Calculate approximate parameters
            original_means <- rowMeans(as.matrix(target_data))
            corrected_means <- rowMeans(as.matrix(corrected_target))
            original_vars <- apply(as.matrix(target_data), 1, var)
            corrected_vars <- apply(as.matrix(corrected_target), 1, var)
            
            # Estimate transformation parameters
            alpha_est <- sqrt(corrected_vars / (original_vars + 1e-8))
            beta_est <- corrected_means - alpha_est * original_means
            
            # Remove infinite/NaN values
            alpha_est[!is.finite(alpha_est)] <- 1.0
            beta_est[!is.finite(beta_est)] <- 0.0
            
            cat("COMBAT_PARAMS:", mean(alpha_est, na.rm=TRUE), mean(beta_est, na.rm=TRUE), "\\n")
            """
            
            # Execute R script
            result = subprocess.run(['R', '--slave', '--vanilla'], 
                                  input=r_script, text=True, capture_output=True)
            
            if result.returncode != 0:
                raise Exception(f"ComBat R script failed: {result.stderr}")
            
            # Load corrected data
            corrected_data = np.loadtxt('temp_combat_result.txt', delimiter='\t')
            
            # Extract parameters from output
            alpha_mean, beta_mean = 1.0, 0.0
            for line in result.stdout.split('\n'):
                if line.startswith('COMBAT_PARAMS:'):
                    parts = line.split()
                    if len(parts) >= 3:
                        alpha_mean = float(parts[1])
                        beta_mean = float(parts[2])
            
            # Clean up temp files
            for temp_file in ['temp_ref_data.txt', 'temp_target_data.txt', 'temp_combat_result.txt']:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
            
            # Create corrected BatchData (ensure correct shape)
            if corrected_data.ndim == 1:
                corrected_data = corrected_data.reshape(-1, 1)
            
            corrected_batch = BatchData(data=corrected_data, gene_indices=target_data.gene_indices)
            metadata = {'alpha_mean': alpha_mean, 'beta_mean': beta_mean}
            
            return corrected_batch, metadata
            
        except Exception as e:
            # Clean up temp files on error
            for temp_file in ['temp_ref_data.txt', 'temp_target_data.txt', 'temp_combat_result.txt']:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
            raise Exception(f"ComBat correction failed: {e}")

    def run_gmm_correction(self, ref_data, target_data):
        """Run GMM correction (simplified implementation)"""
        try:
            from sklearn.mixture import GaussianMixture
            
            # Simple GMM-based correction
            # Fit GMM to reference data
            ref_samples = ref_data.data.T  # Samples x genes
            target_samples = target_data.data.T
            
            # Use subset of genes for speed
            n_genes_subset = min(1000, ref_samples.shape[1])
            gene_indices = np.random.choice(ref_samples.shape[1], n_genes_subset, replace=False)
            
            ref_subset = ref_samples[:, gene_indices]
            target_subset = target_samples[:, gene_indices]
            
            # Fit GMM to reference
            gmm_ref = GaussianMixture(n_components=2, random_state=42)
            gmm_ref.fit(ref_subset)
            
            # Fit GMM to target
            gmm_target = GaussianMixture(n_components=2, random_state=42)
            gmm_target.fit(target_subset)
            
            # Simple correction: match means and variances
            ref_mean = np.mean(ref_samples, axis=0)
            ref_std = np.std(ref_samples, axis=0)
            target_mean = np.mean(target_samples, axis=0)
            target_std = np.std(target_samples, axis=0)
            
            # Standardize target to match reference
            alpha_est = ref_std / (target_std + 1e-8)
            beta_est = ref_mean - alpha_est * target_mean
            
            corrected_samples = alpha_est * target_samples + beta_est
            
            # Create corrected BatchData
            corrected_batch = BatchData(data=corrected_samples.T, gene_indices=target_data.gene_indices)
            metadata = {'alpha_mean': np.mean(alpha_est), 'beta_mean': np.mean(beta_est)}
            
            return corrected_batch, metadata
            
        except Exception as e:
            raise Exception(f"GMM correction failed: {e}")

    def run_posse_correction(self, ref_data, target_data):
        """Run POSSE correction"""
        try:
            posse = POSSE(
                pathway_dict=self._get_default_pathways(ref_data.gene_indices),
                hyperparams=POSSEHyperparameters(tau=25.0, top_k_percent=0.2)
            )
            
            with SuppressOutput():
                corrected_data, metadata = posse.align(ref_data, target_data)
            
            return corrected_data, metadata
            
        except Exception as e:
            raise Exception(f"POSSE correction failed: {e}")

    def run_posse_v6_correction(self, ref_data, target_data):
        """Run POSSE v6.0 correction"""
        try:
            # Convert to v6 BatchData format
            ref_v6 = BatchDataV6(data=ref_data.data, gene_indices=ref_data.gene_indices)
            target_v6 = BatchDataV6(data=target_data.data, gene_indices=target_data.gene_indices)
            
            posse = POSSEv6(
                pathway_dict=self._get_default_pathways(ref_data.gene_indices),
                hyperparams=POSSEv6Hyperparameters(tau=25.0, top_k_percent=0.2)
            )
            
            corrected_data, metadata = posse.align(ref_v6, target_v6)
            
            # Convert back to original BatchData format
            result = BatchData(data=corrected_data.data, gene_indices=corrected_data.gene_indices)
            return result, metadata
            
        except Exception as e:
            raise Exception(f"POSSE v6 correction failed: {e}")

    def silence_test_all_methods(self, dat_lst, gene_names):
        """Test all methods on identical datasets"""
        print("🔬 Comparative Silence Test")
        
        results = {
            'POSSE': {'rmse': None, 'alpha_dev': None, 'beta_dev': None, 'passed': False},
            'POSSEv6': {'rmse': None, 'alpha_dev': None, 'beta_dev': None, 'passed': False},
            'ComBat': {'rmse': None, 'alpha_dev': None, 'beta_dev': None, 'passed': False},
            'GMM': {'rmse': None, 'alpha_dev': None, 'beta_dev': None, 'passed': False}
        }
        
        # Use first study for silence test
        study_names = list(dat_lst.keys())
        test_study = study_names[0]
        print(f"   Testing with study: {test_study}")
        
        original_data = BatchData(data=dat_lst[test_study], gene_indices=gene_names)
        identical_copy = BatchData(data=dat_lst[test_study].copy(), gene_indices=gene_names)
        
        methods = {
            'POSSE': self.run_posse_correction,
            'POSSEv6': self.run_posse_v6_correction,
            'ComBat': self.run_combat_correction,
            'GMM': self.run_gmm_correction
        }
        
        for method_name, method_func in methods.items():
            try:
                corrected_data, metadata = method_func(original_data, identical_copy)
                
                # Calculate metrics
                rmse = np.sqrt(np.mean((corrected_data.data - identical_copy.data)**2))
                alpha_dev = abs(metadata.get('alpha_mean', 1.0) - 1.0)
                beta_dev = abs(metadata.get('beta_mean', 0.0) - 0.0)
                
                # Check if passed (realistic thresholds)
                passed = (rmse < self.config.silence_rmse_threshold and 
                         alpha_dev < self.config.alpha_deviation_threshold and
                         beta_dev < 0.5)
                
                results[method_name] = {
                    'rmse': rmse,
                    'alpha_dev': alpha_dev,
                    'beta_dev': beta_dev,
                    'passed': passed
                }
                
                status = "✅ PASS" if passed else "❌ FAIL"
                print(f"   {method_name}: {status} RMSE={rmse:.3e}, α_dev={alpha_dev:.3f}, β_dev={beta_dev:.3f}")
                
                if not passed:
                    self.method_failures[method_name].append(f"Silence test failed")
                
            except Exception as e:
                print(f"   {method_name}: ❌ CRASH - {e}")
                results[method_name]['passed'] = False
                self.method_failures[method_name].append(f"Silence test crashed: {e}")
        
        return results

    def injection_test_all_methods(self, dat_lst, gene_names):
        """Test all methods on realistic gene-specific synthetic injection"""
        print("🔬 Comparative Injection Test (Gene-Specific Batch Effects)")
        
        results = {
            'POSSE': {},
            'POSSEv6': {},
            'ComBat': {},
            'GMM': {}
        }
        
        # Use first study
        study_names = list(dat_lst.keys())
        test_study = study_names[0]
        original_data = dat_lst[test_study]
        
        print(f"   Testing realistic gene-specific batch effects on: {test_study}")
        
        # Create realistic gene-specific batch effects
        n_genes = original_data.shape[0]
        
        # Global component (dominant signal)
        global_alpha = 1.2
        global_beta = 0.5
        
        # Gene-specific component (realistic heterogeneity)
        np.random.seed(42)  # Reproducible
        gene_alpha = np.random.normal(1.0, 0.2, size=(n_genes, 1))  # ±20% variation
        gene_beta = np.random.normal(0.0, 0.5, size=(n_genes, 1))   # ±0.5 shift variation
        
        # Combine global and gene-specific effects
        final_alpha = global_alpha * gene_alpha
        final_beta = global_beta + gene_beta
        
        # Apply gene-specific distortion
        distorted_data = (original_data * final_alpha) + final_beta
        
        print(f"   Applied gene-specific batch effects:")
        print(f"     Global: α={global_alpha}, β={global_beta}")
        print(f"     Gene-specific α range: [{np.min(final_alpha):.3f}, {np.max(final_alpha):.3f}]")
        print(f"     Gene-specific β range: [{np.min(final_beta):.3f}, {np.max(final_beta):.3f}]")
        
        ref_batch = BatchData(data=original_data, gene_indices=gene_names)
        target_batch = BatchData(data=distorted_data, gene_indices=gene_names)
        
        methods = {
            'POSSE': self.run_posse_correction,
            'POSSEv6': self.run_posse_v6_correction,
            'ComBat': self.run_combat_correction,
            'GMM': self.run_gmm_correction
        }
        
        for method_name, method_func in methods.items():
            try:
                corrected_data, metadata = method_func(ref_batch, target_batch)
                
                # Calculate recovery errors (gene-wise)
                recovery_rmse = np.sqrt(np.mean((corrected_data.data - original_data)**2))
                
                # Calculate parameter recovery (compare to expected inverse transformation)
                estimated_alpha = metadata.get('alpha_mean', 1.0)
                estimated_beta = metadata.get('beta_mean', 0.0)
                
                # For gene-specific effects, we expect approximate recovery of global inverse
                expected_alpha_global = 1.0 / global_alpha  # ≈ 0.833
                expected_beta_global = -global_beta / global_alpha  # ≈ -0.417
                
                alpha_error = abs(estimated_alpha - expected_alpha_global)
                beta_error = abs(estimated_beta - expected_beta_global)
                
                # More lenient thresholds for gene-specific effects
                passed = (alpha_error < 0.5 and beta_error < 0.5 and recovery_rmse < 1.0)
                
                results[method_name] = {
                    'alpha_error': alpha_error,
                    'beta_error': beta_error,
                    'recovery_rmse': recovery_rmse,
                    'expected_alpha': expected_alpha_global,
                    'expected_beta': expected_beta_global,
                    'estimated_alpha': estimated_alpha,
                    'estimated_beta': estimated_beta,
                    'passed': passed
                }
                
                status = "✅ PASS" if passed else "❌ FAIL"
                print(f"   {method_name}: {status}")
                print(f"     Expected: α={expected_alpha_global:.3f}, β={expected_beta_global:.3f}")
                print(f"     Estimated: α={estimated_alpha:.3f}, β={estimated_beta:.3f}")
                print(f"     Errors: α_err={alpha_error:.3f}, β_err={beta_error:.3f}, RMSE={recovery_rmse:.3f}")
                
                if not passed:
                    self.method_failures[method_name].append(f"Gene-specific injection recovery failed")
                
            except Exception as e:
                print(f"   {method_name}: ❌ CRASH - {e}")
                results[method_name] = {'passed': False}
                self.method_failures[method_name].append(f"Gene-specific injection test crashed: {e}")
        
        return results

    def real_data_test_all_methods(self, dat_lst, label_lst, gene_names):
        """Test all methods on real cross-study data"""
        print("🔬 Comparative Real Data Test")
        
        results = {
            'POSSE': {},
            'POSSEv6': {},
            'ComBat': {},
            'GMM': {}
        }
        
        # Test on one study pair
        study_names = list(dat_lst.keys())
        ref_study = study_names[0]  # Africa
        target_study = study_names[1]  # GSE37250_M
        
        print(f"   Testing: {ref_study} → {target_study}")
        
        ref_data = BatchData(data=dat_lst[ref_study], gene_indices=gene_names)
        target_data = BatchData(data=dat_lst[target_study], gene_indices=gene_names)
        
        methods = {
            'POSSE': self.run_posse_correction,
            'POSSEv6': self.run_posse_v6_correction,
            'ComBat': self.run_combat_correction,
            'GMM': self.run_gmm_correction
        }
        
        for method_name, method_func in methods.items():
            try:
                corrected_data, metadata = method_func(ref_data, target_data)
                
                # Calculate basic metrics - handle different metadata key names
                alpha_mean = metadata.get('alpha_final_mean', metadata.get('alpha_mean', 1.0))
                beta_mean = metadata.get('beta_final_mean', metadata.get('beta_mean', 0.0))
                
                # Calculate variance reduction (batch effect removal indicator)
                original_var = np.var(target_data.data, axis=1)
                corrected_var = np.var(corrected_data.data, axis=1)
                var_change = np.mean(corrected_var / (original_var + 1e-8))
                
                results[method_name] = {
                    'alpha_mean': alpha_mean,
                    'beta_mean': beta_mean,
                    'variance_change': var_change,
                    'succeeded': True
                }
                
                print(f"   {method_name}: ✅ SUCCESS α={alpha_mean:.3f}, β={beta_mean:.3f}, var_ratio={var_change:.3f}")
                
            except Exception as e:
                print(f"   {method_name}: ❌ CRASH - {e}")
                results[method_name] = {'succeeded': False}
                self.method_failures[method_name].append(f"Real data test crashed: {e}")
        
        return results

    def topology_stress_test_all_methods(self, dat_lst, label_lst, gene_names):
        """Test all methods on topology preservation (bimodality)"""
        print("🔬 Comparative Topology Stress Test (Bimodality Preservation)")
        
        results = {
            'POSSE': {},
            'POSSEv6': {},
            'ComBat': {},
            'GMM': {}
        }
        
        def bimodality_coefficient(data):
            """Calculate Sarle's bimodality coefficient"""
            if len(data) < 4:
                return 0.0
            
            clean_data = data[np.isfinite(data)]
            if len(clean_data) < 4:
                return 0.0
                
            skewness = stats.skew(clean_data)
            kurtosis = stats.kurtosis(clean_data, fisher=True)
            n = len(clean_data)
            
            bc = (skewness**2 + 1) / (kurtosis + 3 * (n-1)**2 / ((n-2)*(n-3)))
            return bc
        
        # Test genes likely to be bimodal in TB data
        candidate_genes = ['GBP1', 'GBP2', 'IFNG', 'TNF', 'IL1B']
        test_genes = [g for g in candidate_genes if g in gene_names]
        
        if len(test_genes) == 0:
            print("   ⚠️  No candidate bimodal genes found")
            return results
        
        print(f"   Testing bimodality preservation for: {test_genes}")
        
        # Test on first two studies
        study_names = list(dat_lst.keys())
        if len(study_names) < 2:
            print("   ⚠️  Need at least 2 studies for topology test")
            return results
            
        ref_study = study_names[0]
        target_study = study_names[1]
        
        print(f"   Testing: {ref_study} → {target_study}")
        
        ref_data = BatchData(data=dat_lst[ref_study], gene_indices=gene_names)
        target_data = BatchData(data=dat_lst[target_study], gene_indices=gene_names)
        
        methods = {
            'POSSE': self.run_posse_correction,
            'POSSEv6': self.run_posse_v6_correction,
            'ComBat': self.run_combat_correction,
            'GMM': self.run_gmm_correction
        }
        
        for method_name, method_func in methods.items():
            try:
                corrected_data, metadata = method_func(ref_data, target_data)
                
                topology_violations = 0
                gene_results = {}
                
                for gene in test_genes:
                    gene_idx = np.where(gene_names == gene)[0]
                    if len(gene_idx) == 0:
                        continue
                    gene_idx = gene_idx[0]
                    
                    # Get gene expression before and after correction
                    original_expr = target_data.data[gene_idx, :]
                    corrected_expr = corrected_data.data[gene_idx, :]
                    
                    # Calculate bimodality coefficients
                    original_bc = bimodality_coefficient(original_expr)
                    corrected_bc = bimodality_coefficient(corrected_expr)
                    
                    # Calculate topology preservation
                    bc_change = abs(corrected_bc - original_bc)
                    topology_preserved = bc_change < self.config.topology_preservation_threshold
                    
                    gene_results[gene] = {
                        'original_bc': original_bc,
                        'corrected_bc': corrected_bc,
                        'bc_change': bc_change,
                        'preserved': topology_preserved
                    }
                    
                    if not topology_preserved:
                        topology_violations += 1
                
                # Overall method assessment
                total_genes = len(gene_results)
                preservation_rate = (total_genes - topology_violations) / max(1, total_genes)
                passed = topology_violations == 0
                
                results[method_name] = {
                    'gene_results': gene_results,
                    'topology_violations': topology_violations,
                    'total_genes': total_genes,
                    'preservation_rate': preservation_rate,
                    'passed': passed
                }
                
                status = "✅ PASS" if passed else "❌ FAIL"
                print(f"   {method_name}: {status}")
                print(f"     Topology preserved: {total_genes - topology_violations}/{total_genes} genes ({preservation_rate:.1%})")
                
                if not passed:
                    self.method_failures[method_name].append(f"Topology violations: {topology_violations}/{total_genes} genes")
                
            except Exception as e:
                print(f"   {method_name}: ❌ CRASH - {e}")
                results[method_name] = {'passed': False}
                self.method_failures[method_name].append(f"Topology test crashed: {e}")
        
        return results

    def _get_default_pathways(self, gene_names):
        """Get default pathway dictionary filtered for available genes"""
        default_pathways = {
            'interferon': ['IFNG', 'IRF1', 'IRF7', 'STAT1', 'GBP1', 'GBP2', 'OAS1', 'MX1'],
            'inflammatory': ['TNF', 'IL1B', 'IL6', 'NFKB1', 'RELA', 'PTGS2', 'ICAM1'],
            'housekeeping': ['ACTB', 'GAPDH', 'B2M', 'HPRT1', 'TBP', 'YWHAZ']
        }
        
        filtered_pathways = {}
        for pathway_name, pathway_genes in default_pathways.items():
            present_genes = [g for g in pathway_genes if g in gene_names]
            if len(present_genes) >= 3:
                filtered_pathways[pathway_name] = present_genes
        
        return filtered_pathways

    def run_comparative_validation(self, dat_lst, label_lst, gene_names):
        """Run complete comparative validation"""
        print("🚀 Comparative Validation Suite - POSSE vs ComBat vs GMM")
        print("="*70)
        
        # Clear previous failures
        self.method_failures = {'POSSE': [], 'POSSEv6': [], 'ComBat': [], 'GMM': []}
        
        # Run tests
        silence_results = self.silence_test_all_methods(dat_lst, gene_names)
        print()
        
        injection_results = self.injection_test_all_methods(dat_lst, gene_names)
        print()
        
        real_data_results = self.real_data_test_all_methods(dat_lst, label_lst, gene_names)
        print()
        
        topology_results = self.topology_stress_test_all_methods(dat_lst, label_lst, gene_names)
        print()
        
        # Store results
        self.results = {
            'silence_test': silence_results,
            'injection_test': injection_results,
            'real_data_test': real_data_results,
            'topology_test': topology_results
        }
        
        # Generate comparative report
        self.generate_comparative_report()
        
        return self.results

    def generate_comparative_report(self):
        """Generate comparative validation report"""
        print("📊 COMPARATIVE VALIDATION REPORT")
        print("="*70)
        
        methods = ['POSSE', 'POSSEv6', 'ComBat', 'GMM']
        
        # Count failures per method
        failure_counts = {method: len(self.method_failures[method]) for method in methods}
        
        print("🏆 METHOD RANKINGS (by failure count):")
        sorted_methods = sorted(methods, key=lambda m: failure_counts[m])
        
        for i, method in enumerate(sorted_methods, 1):
            failures = failure_counts[method]
            status = "🥇" if i == 1 else "🥈" if i == 2 else "🥉"
            print(f"   {status} {method}: {failures} failures")
        
        print()
        print("📋 DETAILED RESULTS:")
        
        # Silence test summary
        if 'silence_test' in self.results:
            print("   Silence Test (Identical Datasets):")
            for method in methods:
                result = self.results['silence_test'][method]
                if result.get('rmse') is not None:
                    status = "✅" if result['passed'] else "❌"
                    print(f"     {method}: {status} RMSE={result['rmse']:.3e}")
                else:
                    print(f"     {method}: ❌ CRASHED")
        
        # Injection test summary
        if 'injection_test' in self.results:
            print("   Injection Test (Ground Truth Recovery):")
            for method in methods:
                result = self.results['injection_test'][method]
                if result.get('alpha_error') is not None:
                    status = "✅" if result['passed'] else "❌"
                    print(f"     {method}: {status} α_err={result['alpha_error']:.3f}")
                else:
                    print(f"     {method}: ❌ CRASHED")
        
        # Real data test summary
        if 'real_data_test' in self.results:
            print("   Real Data Test (Cross-Study Correction):")
            for method in methods:
                result = self.results['real_data_test'][method]
                if result.get('succeeded'):
                    print(f"     {method}: ✅ α={result['alpha_mean']:.3f}")
                else:
                    print(f"     {method}: ❌ CRASHED")
        
        # Topology test summary
        if 'topology_test' in self.results:
            print("   Topology Test (Bimodality Preservation):")
            for method in methods:
                result = self.results['topology_test'][method]
                if result.get('passed') is not None:
                    status = "✅" if result['passed'] else "❌"
                    rate = result.get('preservation_rate', 0)
                    print(f"     {method}: {status} {rate:.1%} preserved")
                else:
                    print(f"     {method}: ❌ CRASHED")
        
        print()
        print("🎯 COMPARATIVE INSIGHTS:")
        
        best_method = sorted_methods[0]
        worst_method = sorted_methods[-1]
        
        print(f"   Best performing method: {best_method} ({failure_counts[best_method]} failures)")
        print(f"   Worst performing method: {worst_method} ({failure_counts[worst_method]} failures)")
        
        # Method-specific insights
        if failure_counts['POSSE'] > failure_counts['ComBat']:
            print("   POSSE underperforms ComBat - needs fundamental fixes")
        elif failure_counts['POSSE'] < failure_counts['ComBat']:
            print("   POSSE outperforms ComBat - pathway approach shows promise")
        else:
            print("   POSSE and ComBat show similar performance")
        
        if failure_counts['GMM'] < min(failure_counts['POSSE'], failure_counts['ComBat']):
            print("   GMM shows superior robustness for this dataset")
        
        print()
        print("📈 RECOMMENDATIONS:")
        if failure_counts[best_method] == 0:
            print(f"   {best_method} passes all tests - ready for production")
        elif failure_counts[best_method] <= 2:
            print(f"   {best_method} shows good performance - minor fixes needed")
        else:
            print("   All methods show significant issues - batch correction is inherently difficult")
            print("   Focus on the method with most promise for improvement")

def main():
    """Main function"""
    print("Comparative Validation Suite")
    print("Tests POSSE vs ComBat vs GMM with realistic expectations")
    
    # Placeholder for real data
    dat_lst = {}
    label_lst = {}
    gene_names = np.array([])
    
    print("Framework ready. Load real data to execute comparative validation.")

if __name__ == "__main__":
    main()