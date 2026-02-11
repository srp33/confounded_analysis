#!/usr/bin/env python3
"""
POSSE Validation Suite v2.0
Rigorous method validation following the diagnostic critique recommendations
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

# Import existing modules
from pace import BatchData
from posse import POSSE, POSSEHyperparameters

@dataclass
class ValidationConfig:
    """Configuration for validation experiments"""
    output_dir: str = "posse_validation"
    save_plots: bool = True
    verbose: bool = True
    random_state: int = 42
    # Failure thresholds
    trust_threshold: float = 0.5
    alpha_deviation_threshold: float = 0.2
    silence_rmse_threshold: float = 1e-9
    topology_preservation_threshold: float = 0.1

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

class POSSEValidationSuite:
    """Rigorous validation suite for POSSE method validation"""
    
    def __init__(self, config: ValidationConfig = None):
        self.config = config or ValidationConfig()
        self.results = {}
        self.failures = []
        
    def experiment_0_sanity_check(self, dat_lst, label_lst, gene_names):
        """
        Experiment 0: Sanity Check (Silence Test)
        
        Tests if POSSE correctly handles identical datasets:
        - Load Study A, create Study A' (exact copy)
        - Run POSSE alignment
        - Assert α=1.0, β=0.0, RMSE<1e-9
        """
        print("🔬 Experiment 0: Sanity Check (Silence Test)")
        
        results = {
            'silence_tests': {},
            'failures': [],
            'rmse_values': {},
            'parameter_deviations': {}
        }
        
        # Test with first available study
        study_names = list(dat_lst.keys())
        if len(study_names) == 0:
            print("   ❌ No studies available for silence test")
            return results
            
        test_study = study_names[0]
        print(f"   Testing silence with study: {test_study}")
        
        # Create identical copy
        original_data = BatchData(data=dat_lst[test_study], gene_indices=gene_names)
        identical_copy = BatchData(data=dat_lst[test_study].copy(), gene_indices=gene_names)
        
        # Test POSSE on identical data
        posse = POSSE(
            pathway_dict=self._get_default_pathways(gene_names),
            hyperparams=POSSEHyperparameters(tau=25.0, top_k_percent=0.2)
        )
        
        try:
            with SuppressOutput():
                corrected_data, metadata = posse.align(original_data, identical_copy)
            
            # Calculate deviations
            alpha_mean = metadata.get('alpha_mean', 1.0)
            beta_mean = metadata.get('beta_mean', 0.0)
            
            # Calculate RMSE between original and "corrected" identical data
            rmse = np.sqrt(np.mean((corrected_data.data - identical_copy.data)**2))
            
            # Check silence criteria
            alpha_deviation = abs(alpha_mean - 1.0)
            beta_deviation = abs(beta_mean - 0.0)
            
            results['silence_tests'][test_study] = {
                'alpha_mean': alpha_mean,
                'beta_mean': beta_mean,
                'rmse': rmse,
                'alpha_deviation': alpha_deviation,
                'beta_deviation': beta_deviation
            }
            
            # Flag failures
            silence_pass = True
            if rmse > self.config.silence_rmse_threshold:
                failure = f"RMSE too high: {rmse:.2e} > {self.config.silence_rmse_threshold:.2e}"
                results['failures'].append(failure)
                self.failures.append(f"Silence Test: {failure}")
                silence_pass = False
                
            if alpha_deviation > self.config.alpha_deviation_threshold:
                failure = f"Alpha deviation: {alpha_deviation:.3f} > {self.config.alpha_deviation_threshold:.3f}"
                results['failures'].append(failure)
                self.failures.append(f"Silence Test: {failure}")
                silence_pass = False
                
            if beta_deviation > 0.1:  # Beta should be very close to 0
                failure = f"Beta deviation: {beta_deviation:.3f} > 0.1"
                results['failures'].append(failure)
                self.failures.append(f"Silence Test: {failure}")
                silence_pass = False
            
            status = "✅ PASS" if silence_pass else "❌ FAIL"
            print(f"   {status} α={alpha_mean:.6f}, β={beta_mean:.6f}, RMSE={rmse:.2e}")
            
        except Exception as e:
            failure = f"Silence test crashed: {e}"
            results['failures'].append(failure)
            self.failures.append(f"Silence Test: {failure}")
            print(f"   ❌ FAIL: {failure}")
        
        self.results['sanity_check'] = results
        return results

    def experiment_1_synthetic_injection(self, dat_lst, label_lst, gene_names):
        """
        Experiment 1: Synthetic Injection (Ground Truth Recovery)
        
        Tests POSSE's ability to recover known distortions:
        - Take Study A, create Study B = Study A + known shift/scale
        - Run POSSE to correct B → A
        - Measure recovery accuracy
        """
        print("🔬 Experiment 1: Synthetic Injection (Ground Truth Recovery)")
        
        results = {
            'injection_tests': {},
            'recovery_errors': {},
            'failures': []
        }
        
        study_names = list(dat_lst.keys())
        if len(study_names) == 0:
            print("   ❌ No studies available for injection test")
            return results
            
        test_study = study_names[0]
        print(f"   Testing injection recovery with study: {test_study}")
        
        # Define known distortions to inject
        distortions = [
            {'name': 'shift_+2.0', 'alpha': 1.0, 'beta': 2.0},
            {'name': 'scale_0.5', 'alpha': 0.5, 'beta': 0.0},
            {'name': 'shift_scale', 'alpha': 1.5, 'beta': -1.0}
        ]
        
        original_data = dat_lst[test_study]
        
        for distortion in distortions:
            print(f"   Testing distortion: {distortion['name']}")
            
            # Inject known distortion: Y = α*X + β
            true_alpha = distortion['alpha']
            true_beta = distortion['beta']
            
            distorted_data = true_alpha * original_data + true_beta
            
            # Create BatchData objects
            ref_batch = BatchData(data=original_data, gene_indices=gene_names)
            target_batch = BatchData(data=distorted_data, gene_indices=gene_names)
            
            # Run POSSE to recover original from distorted
            posse = POSSE(
                pathway_dict=self._get_default_pathways(gene_names),
                hyperparams=POSSEHyperparameters(tau=25.0, top_k_percent=0.2)
            )
            
            try:
                with SuppressOutput():
                    corrected_data, metadata = posse.align(ref_batch, target_batch)
                
                # Extract POSSE's estimated parameters
                estimated_alpha = metadata.get('alpha_mean', 1.0)
                estimated_beta = metadata.get('beta_mean', 0.0)
                
                # Calculate recovery errors
                alpha_error = abs(estimated_alpha - (1.0/true_alpha))  # POSSE should find inverse
                beta_error = abs(estimated_beta - (-true_beta/true_alpha))
                
                # Calculate data recovery error
                recovery_rmse = np.sqrt(np.mean((corrected_data.data - original_data)**2))
                
                results['injection_tests'][distortion['name']] = {
                    'true_alpha': true_alpha,
                    'true_beta': true_beta,
                    'estimated_alpha': estimated_alpha,
                    'estimated_beta': estimated_beta,
                    'alpha_error': alpha_error,
                    'beta_error': beta_error,
                    'recovery_rmse': recovery_rmse
                }
                
                # Check recovery quality
                recovery_pass = True
                if alpha_error > 0.1:
                    failure = f"{distortion['name']}: Alpha recovery error {alpha_error:.3f} > 0.1"
                    results['failures'].append(failure)
                    self.failures.append(f"Injection Test: {failure}")
                    recovery_pass = False
                    
                if beta_error > 0.1:
                    failure = f"{distortion['name']}: Beta recovery error {beta_error:.3f} > 0.1"
                    results['failures'].append(failure)
                    self.failures.append(f"Injection Test: {failure}")
                    recovery_pass = False
                
                status = "✅ PASS" if recovery_pass else "❌ FAIL"
                print(f"     {status} α_err={alpha_error:.3f}, β_err={beta_error:.3f}, RMSE={recovery_rmse:.3f}")
                
            except Exception as e:
                failure = f"{distortion['name']}: Injection test crashed: {e}"
                results['failures'].append(failure)
                self.failures.append(f"Injection Test: {failure}")
                print(f"     ❌ FAIL: {failure}")
        
        self.results['synthetic_injection'] = results
        return results

    def experiment_2_frankenstein_monitor(self, dat_lst, label_lst, gene_names):
        """
        Experiment 2: Frankenstein Monitor (Trust vs Deviation Boundary)
        
        Tests POSSE's trust system with explicit failure boundaries:
        - Calculate trust scores and parameter deviations
        - Flag violations: if trust < threshold and |α-1| > threshold
        - Generate Trust vs Deviation diagnostic plots
        """
        print("🔬 Experiment 2: Frankenstein Monitor (Trust Boundary Detection)")
        
        results = {
            'trust_violations': [],
            'trust_data': {},
            'boundary_failures': [],
            'diagnostic_plots': {}
        }
        
        # Create diagnostic POSSE that exposes trust calculations
        class FrankensteinPOSSE(POSSE):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.trust_diagnostics = {}
            
            def pathway_execution(self, X_prime, Y_prime, pathway_indices, C_null, global_gene_idxs, alpha_prior_vec, beta_prior_vec, tau_override=None, top_k_override=None):
                omega, alpha_est, beta_est, k_raw, metrics = super().pathway_execution(
                    X_prime, Y_prime, pathway_indices, C_null, global_gene_idxs, 
                    alpha_prior_vec, beta_prior_vec, tau_override, top_k_override
                )
                
                # Store detailed trust diagnostics
                pathway_name = f"pathway_{len(self.trust_diagnostics)}"
                self.trust_diagnostics[pathway_name] = {
                    'omega': omega,
                    'alpha_estimates': alpha_est,
                    'beta_estimates': beta_est,
                    'trust_scores': metrics.get('avg_trust', 0.0),
                    'correlations': metrics.get('avg_correlation', 0.0),
                    'gene_indices': global_gene_idxs
                }
                
                return omega, alpha_est, beta_est, k_raw, metrics
        
        # Test trust system on real study pairs
        study_names = list(dat_lst.keys())
        
        for i, ref_study in enumerate(study_names[:2]):  # Test first 2 studies
            for target_study in study_names[i+1:i+2]:
                print(f"   Analyzing trust boundary: {ref_study} → {target_study}")
                
                ref_data = BatchData(data=dat_lst[ref_study], gene_indices=gene_names)
                target_data = BatchData(data=dat_lst[target_study], gene_indices=gene_names)
                
                frankenstein_posse = FrankensteinPOSSE(
                    pathway_dict=self._get_default_pathways(gene_names),
                    hyperparams=POSSEHyperparameters(tau=25.0, top_k_percent=0.2)
                )
                
                try:
                    with SuppressOutput():
                        corrected_data, metadata = frankenstein_posse.align(ref_data, target_data)
                    
                    # Analyze trust vs deviation for each pathway
                    trust_violations = []
                    trust_plot_data = {'trust': [], 'alpha_deviation': [], 'pathway': []}
                    
                    for pathway_name, diagnostics in frankenstein_posse.trust_diagnostics.items():
                        trust_score = diagnostics['trust_scores']
                        alpha_estimates = diagnostics['alpha_estimates']
                        
                        # Calculate deviation from identity (α=1.0)
                        alpha_deviations = np.abs(alpha_estimates - 1.0)
                        mean_alpha_deviation = np.mean(alpha_deviations)
                        
                        # Store for plotting
                        trust_plot_data['trust'].append(trust_score)
                        trust_plot_data['alpha_deviation'].append(mean_alpha_deviation)
                        trust_plot_data['pathway'].append(pathway_name)
                        
                        # Check Frankenstein boundary violation
                        if (trust_score < self.config.trust_threshold and 
                            mean_alpha_deviation > self.config.alpha_deviation_threshold):
                            
                            violation = {
                                'pathway': pathway_name,
                                'trust_score': trust_score,
                                'alpha_deviation': mean_alpha_deviation,
                                'study_pair': f"{ref_study}→{target_study}"
                            }
                            trust_violations.append(violation)
                            
                            failure_msg = (f"Trust violation in {pathway_name}: "
                                         f"trust={trust_score:.3f} < {self.config.trust_threshold}, "
                                         f"α_dev={mean_alpha_deviation:.3f} > {self.config.alpha_deviation_threshold}")
                            
                            results['boundary_failures'].append(failure_msg)
                            self.failures.append(f"Frankenstein Monitor: {failure_msg}")
                            print(f"     ❌ VIOLATION: {failure_msg}")
                    
                    results['trust_violations'].extend(trust_violations)
                    results['trust_data'][f"{ref_study}_{target_study}"] = trust_plot_data
                    
                    if len(trust_violations) == 0:
                        print(f"     ✅ No trust boundary violations detected")
                    else:
                        print(f"     ❌ {len(trust_violations)} trust violations found")
                        
                except Exception as e:
                    failure = f"Frankenstein monitor crashed on {ref_study}→{target_study}: {e}"
                    results['boundary_failures'].append(failure)
                    self.failures.append(f"Frankenstein Monitor: {failure}")
                    print(f"     ❌ FAIL: {failure}")
        
        self.results['frankenstein_monitor'] = results
        return results

    def experiment_3_topology_stress(self, dat_lst, label_lst, gene_names):
        """
        Experiment 3: Topology Stress Test (Bimodality Preservation)
        
        Tests if POSSE preserves distributional topology:
        - Identify bimodal genes (like GBP1 in TB data)
        - Measure bimodality coefficient before/after correction
        - Assert topology preservation within threshold
        """
        print("🔬 Experiment 3: Topology Stress Test (Bimodality Preservation)")
        
        results = {
            'bimodality_tests': {},
            'topology_failures': [],
            'gene_topology_scores': {}
        }
        
        def bimodality_coefficient(data):
            """Calculate bimodality coefficient (Sarle's bimodality coefficient)"""
            if len(data) < 4:
                return 0.0
            
            # Remove NaN/inf values
            clean_data = data[np.isfinite(data)]
            if len(clean_data) < 4:
                return 0.0
                
            skewness = stats.skew(clean_data)
            kurtosis = stats.kurtosis(clean_data, fisher=True)  # Excess kurtosis
            n = len(clean_data)
            
            # Sarle's bimodality coefficient
            bc = (skewness**2 + 1) / (kurtosis + 3 * (n-1)**2 / ((n-2)*(n-3)))
            return bc
        
        # Test genes likely to be bimodal in TB data
        candidate_genes = ['GBP1', 'GBP2', 'IFNG', 'TNF', 'IL1B', 'IRF1', 'STAT1']
        test_genes = [g for g in candidate_genes if g in gene_names]
        
        if len(test_genes) == 0:
            print("   ❌ No candidate bimodal genes found in dataset")
            return results
        
        print(f"   Testing topology preservation for genes: {test_genes}")
        
        # Test on study pairs
        study_names = list(dat_lst.keys())
        
        for i, ref_study in enumerate(study_names[:2]):
            for target_study in study_names[i+1:i+2]:
                print(f"   Testing topology: {ref_study} → {target_study}")
                
                ref_data = BatchData(data=dat_lst[ref_study], gene_indices=gene_names)
                target_data = BatchData(data=dat_lst[target_study], gene_indices=gene_names)
                
                # Run POSSE correction
                posse = POSSE(
                    pathway_dict=self._get_default_pathways(gene_names),
                    hyperparams=POSSEHyperparameters(tau=25.0, top_k_percent=0.2)
                )
                
                try:
                    with SuppressOutput():
                        corrected_data, metadata = posse.align(ref_data, target_data)
                    
                    topology_violations = []
                    
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
                        
                        results['gene_topology_scores'][f"{gene}_{ref_study}_{target_study}"] = {
                            'gene': gene,
                            'original_bimodality': original_bc,
                            'corrected_bimodality': corrected_bc,
                            'bimodality_change': bc_change,
                            'topology_preserved': topology_preserved
                        }
                        
                        if not topology_preserved:
                            violation = (f"Topology violation for {gene}: "
                                       f"bimodality change {bc_change:.3f} > {self.config.topology_preservation_threshold}")
                            topology_violations.append(violation)
                            results['topology_failures'].append(violation)
                            self.failures.append(f"Topology Stress: {violation}")
                            print(f"     ❌ {violation}")
                        else:
                            print(f"     ✅ {gene}: BC change {bc_change:.3f} (preserved)")
                    
                    results['bimodality_tests'][f"{ref_study}_{target_study}"] = {
                        'violations': topology_violations,
                        'genes_tested': test_genes,
                        'violations_count': len(topology_violations)
                    }
                    
                except Exception as e:
                    failure = f"Topology test crashed on {ref_study}→{target_study}: {e}"
                    results['topology_failures'].append(failure)
                    self.failures.append(f"Topology Stress: {failure}")
                    print(f"     ❌ FAIL: {failure}")
        
        self.results['topology_stress'] = results
        return results

    def experiment_4_compositional_imbalance(self, dat_lst, label_lst, gene_names):
        """
        Experiment 4: Compositional Imbalance Test
        
        Tests if POSSE preserves biological composition differences:
        - Measure cluster alignment before/after correction
        - Check if population-specific signatures are preserved
        - Flag if correction crushes meaningful biological differences
        """
        print("🔬 Experiment 4: Compositional Imbalance Test")
        
        results = {
            'composition_tests': {},
            'cluster_alignment_errors': {},
            'composition_failures': []
        }
        
        # Test composition preservation across studies with different TB severity profiles
        study_names = list(dat_lst.keys())
        
        for i, ref_study in enumerate(study_names[:2]):
            for target_study in study_names[i+1:i+2]:
                print(f"   Testing composition preservation: {ref_study} → {target_study}")
                
                ref_data = BatchData(data=dat_lst[ref_study], gene_indices=gene_names)
                target_data = BatchData(data=dat_lst[target_study], gene_indices=gene_names)
                ref_labels = np.array(label_lst[ref_study])
                target_labels = np.array(label_lst[target_study])
                
                # Calculate original composition (TB vs control ratios)
                ref_tb_ratio = np.mean(ref_labels == 1) if len(ref_labels) > 0 else 0.5
                target_tb_ratio = np.mean(target_labels == 1) if len(target_labels) > 0 else 0.5
                
                print(f"     Original TB ratios: {ref_study}={ref_tb_ratio:.2f}, {target_study}={target_tb_ratio:.2f}")
                
                # Run POSSE correction
                posse = POSSE(
                    pathway_dict=self._get_default_pathways(gene_names),
                    hyperparams=POSSEHyperparameters(tau=25.0, top_k_percent=0.2)
                )
                
                try:
                    with SuppressOutput():
                        corrected_data, metadata = posse.align(ref_data, target_data)
                    
                    # Measure cluster separation before/after correction
                    # Use key TB biomarkers for cluster analysis
                    tb_biomarkers = ['IFNG', 'GBP1', 'TNF', 'IL1B']
                    biomarker_indices = []
                    for gene in tb_biomarkers:
                        idx = np.where(gene_names == gene)[0]
                        if len(idx) > 0:
                            biomarker_indices.append(idx[0])
                    
                    if len(biomarker_indices) >= 2:
                        # Original cluster separation
                        original_biomarkers = target_data.data[biomarker_indices, :]
                        tb_samples = target_labels == 1
                        control_samples = target_labels == 0
                        
                        if np.sum(tb_samples) > 0 and np.sum(control_samples) > 0:
                            original_tb_mean = np.mean(original_biomarkers[:, tb_samples], axis=1)
                            original_control_mean = np.mean(original_biomarkers[:, control_samples], axis=1)
                            original_separation = np.linalg.norm(original_tb_mean - original_control_mean)
                            
                            # Corrected cluster separation
                            corrected_biomarkers = corrected_data.data[biomarker_indices, :]
                            corrected_tb_mean = np.mean(corrected_biomarkers[:, tb_samples], axis=1)
                            corrected_control_mean = np.mean(corrected_biomarkers[:, control_samples], axis=1)
                            corrected_separation = np.linalg.norm(corrected_tb_mean - corrected_control_mean)
                            
                            # Calculate cluster alignment error
                            separation_ratio = corrected_separation / (original_separation + 1e-8)
                            cluster_alignment_error = abs(1.0 - separation_ratio)
                            
                            results['cluster_alignment_errors'][f"{ref_study}_{target_study}"] = {
                                'original_separation': original_separation,
                                'corrected_separation': corrected_separation,
                                'separation_ratio': separation_ratio,
                                'alignment_error': cluster_alignment_error
                            }
                            
                            # Check if biological signal was crushed
                            if cluster_alignment_error > 0.5:  # More than 50% change in separation
                                failure = (f"Cluster alignment error {cluster_alignment_error:.3f} > 0.5 "
                                         f"(separation ratio: {separation_ratio:.3f})")
                                results['composition_failures'].append(failure)
                                self.failures.append(f"Compositional Imbalance: {failure}")
                                print(f"     ❌ {failure}")
                            else:
                                print(f"     ✅ Cluster separation preserved (ratio: {separation_ratio:.3f})")
                        else:
                            print(f"     ⚠️  Insufficient TB/control samples for cluster analysis")
                    else:
                        print(f"     ⚠️  Insufficient biomarker genes for cluster analysis")
                    
                    results['composition_tests'][f"{ref_study}_{target_study}"] = {
                        'ref_tb_ratio': ref_tb_ratio,
                        'target_tb_ratio': target_tb_ratio,
                        'biomarkers_tested': len(biomarker_indices)
                    }
                    
                except Exception as e:
                    failure = f"Composition test crashed on {ref_study}→{target_study}: {e}"
                    results['composition_failures'].append(failure)
                    self.failures.append(f"Compositional Imbalance: {failure}")
                    print(f"     ❌ FAIL: {failure}")
        
        self.results['compositional_imbalance'] = results
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

    def run_validation_suite(self, dat_lst, label_lst, gene_names):
        """Run complete validation suite"""
        print("🚀 POSSE Validation Suite v2.0 - Method Validation")
        print("="*60)
        
        # Clear previous failures
        self.failures = []
        
        # Run validation experiments
        exp0_results = self.experiment_0_sanity_check(dat_lst, label_lst, gene_names)
        print()
        
        exp1_results = self.experiment_1_synthetic_injection(dat_lst, label_lst, gene_names)
        print()
        
        exp2_results = self.experiment_2_frankenstein_monitor(dat_lst, label_lst, gene_names)
        print()
        
        exp3_results = self.experiment_3_topology_stress(dat_lst, label_lst, gene_names)
        print()
        
        exp4_results = self.experiment_4_compositional_imbalance(dat_lst, label_lst, gene_names)
        print()
        
        # Generate validation report
        self.generate_validation_report()
        
        return self.results

    def generate_validation_report(self):
        """Generate comprehensive validation report"""
        print("📊 VALIDATION REPORT")
        print("="*60)
        
        total_failures = len(self.failures)
        
        if total_failures == 0:
            print("🎉 ALL TESTS PASSED - POSSE validation successful!")
        else:
            print(f"❌ {total_failures} VALIDATION FAILURES DETECTED:")
            for i, failure in enumerate(self.failures, 1):
                print(f"   {i}. {failure}")
        
        print()
        print("📋 VALIDATION SUMMARY:")
        
        # Sanity check summary
        if 'sanity_check' in self.results:
            sc = self.results['sanity_check']
            sc_failures = len(sc.get('failures', []))
            print(f"   Sanity Check (Silence Test): {'✅ PASS' if sc_failures == 0 else f'❌ {sc_failures} failures'}")
        
        # Injection test summary
        if 'synthetic_injection' in self.results:
            si = self.results['synthetic_injection']
            si_failures = len(si.get('failures', []))
            print(f"   Synthetic Injection: {'✅ PASS' if si_failures == 0 else f'❌ {si_failures} failures'}")
        
        # Frankenstein monitor summary
        if 'frankenstein_monitor' in self.results:
            fm = self.results['frankenstein_monitor']
            fm_failures = len(fm.get('boundary_failures', []))
            print(f"   Frankenstein Monitor: {'✅ PASS' if fm_failures == 0 else f'❌ {fm_failures} violations'}")
        
        # Topology stress summary
        if 'topology_stress' in self.results:
            ts = self.results['topology_stress']
            ts_failures = len(ts.get('topology_failures', []))
            print(f"   Topology Stress Test: {'✅ PASS' if ts_failures == 0 else f'❌ {ts_failures} violations'}")
        
        # Compositional imbalance summary
        if 'compositional_imbalance' in self.results:
            ci = self.results['compositional_imbalance']
            ci_failures = len(ci.get('composition_failures', []))
            print(f"   Compositional Imbalance: {'✅ PASS' if ci_failures == 0 else f'❌ {ci_failures} violations'}")
        
        print()
        print("🎯 VALIDATION VERDICT:")
        if total_failures == 0:
            print("   POSSE method validation: ✅ PASSED")
            print("   Method is ready for production use")
        else:
            print("   POSSE method validation: ❌ FAILED")
            print("   Method requires fixes before production use")
            print("   Focus areas for POSSE v6.0 development:")
            
            # Categorize failures for actionable insights
            silence_failures = [f for f in self.failures if 'Silence Test' in f]
            injection_failures = [f for f in self.failures if 'Injection Test' in f]
            trust_failures = [f for f in self.failures if 'Frankenstein Monitor' in f]
            topology_failures = [f for f in self.failures if 'Topology Stress' in f]
            composition_failures = [f for f in self.failures if 'Compositional Imbalance' in f]
            
            if silence_failures:
                print("     - Fix basic parameter estimation (silence test failures)")
            if injection_failures:
                print("     - Improve ground truth recovery accuracy")
            if trust_failures:
                print("     - Recalibrate trust system thresholds and logic")
            if topology_failures:
                print("     - Preserve distributional topology (bimodality)")
            if composition_failures:
                print("     - Prevent biological signal crushing")

def main():
    """Main function to run validation suite"""
    
    print("POSSE Validation Suite v2.0")
    print("Rigorous method validation following diagnostic critique")
    print("Note: This framework needs real TB data to run")
    print("To use: load TB_real_data.RData and populate dat_lst, label_lst, gene_names")
    
    # Initialize validation suite
    config = ValidationConfig(
        output_dir="posse_validation_results",
        save_plots=True,
        verbose=True
    )
    
    validation_suite = POSSEValidationSuite(config)
    
    # Placeholder for real data
    dat_lst = {}
    label_lst = {}
    gene_names = np.array([])
    
    print("Validation framework ready. Load real data to execute validation.")

if __name__ == "__main__":
    main()