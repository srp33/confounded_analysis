#!/usr/bin/env python3
"""
POSSE Diagnostic Experiments Suite
Comprehensive analysis to understand POSSE's cross-study TB prediction failures
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
class ExperimentConfig:
    """Configuration for diagnostic experiments"""
    output_dir: str = "posse_diagnostics"
    save_plots: bool = True
    verbose: bool = True
    random_state: int = 42

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

class POSSEDiagnosticSuite:
    """Comprehensive diagnostic suite for POSSE analysis"""
    
    def __init__(self, config: ExperimentConfig = None):
        self.config = config or ExperimentConfig()
        self.results = {}
        
    def load_tb_data(self):
        """Load TB data for analysis"""
        # Implementation will load real TB data
        pass
    def experiment_1_variance_decomposition(self, dat_lst, label_lst, gene_names):
        """
        Experiment 1: Biological vs Technical Variance Decomposition
        
        Decomposes cross-study variance into:
        1. Technical factors (library size, platform, processing)
        2. Population factors (genetics, demographics) 
        3. Disease factors (TB severity, progression stage)
        """
        print("🔬 Experiment 1: Variance Decomposition Analysis")
        
        results = {
            'within_study_var': {},
            'between_study_var': {},
            'technical_genes': [],
            'biological_genes': [],
            'population_genes': []
        }
        
        # Combine all data
        all_data = []
        all_labels = []
        all_studies = []
        
        for study, data in dat_lst.items():
            all_data.append(data.T)  # Transpose to samples x genes
            all_labels.extend(label_lst[study])
            all_studies.extend([study] * len(label_lst[study]))
        
        combined_data = np.vstack(all_data)
        combined_labels = np.array(all_labels)
        combined_studies = np.array(all_studies)
        
        print(f"   Combined data shape: {combined_data.shape}")
        print(f"   Studies: {np.unique(combined_studies)}")
        
        # Calculate variance components for each gene
        n_genes = combined_data.shape[1]
        
        for g in range(min(1000, n_genes)):  # Analyze subset for speed
            gene_expr = combined_data[:, g]
            
            # Within-study variance (technical + biological noise)
            within_var = 0
            n_within = 0
            
            for study in np.unique(combined_studies):
                study_mask = combined_studies == study
                study_expr = gene_expr[study_mask]
                if len(study_expr) > 1:
                    within_var += np.var(study_expr, ddof=1) * (len(study_expr) - 1)
                    n_within += len(study_expr) - 1
            
            within_var = within_var / max(1, n_within)
            
            # Between-study variance (batch effects + population differences)
            study_means = []
            for study in np.unique(combined_studies):
                study_mask = combined_studies == study
                study_expr = gene_expr[study_mask]
                if len(study_expr) > 0:
                    study_means.append(np.mean(study_expr))
            
            between_var = np.var(study_means, ddof=1) if len(study_means) > 1 else 0
            
            # Store results
            gene_name = gene_names[g] if g < len(gene_names) else f"Gene_{g}"
            results['within_study_var'][gene_name] = within_var
            results['between_study_var'][gene_name] = between_var
            
            # Classify genes based on variance patterns
            if between_var > 2 * within_var:
                results['biological_genes'].append(gene_name)
            elif within_var > 2 * between_var:
                results['technical_genes'].append(gene_name)
        
        # Summary statistics
        within_vars = list(results['within_study_var'].values())
        between_vars = list(results['between_study_var'].values())
        
        print(f"   Within-study variance: {np.mean(within_vars):.3f} ± {np.std(within_vars):.3f}")
        print(f"   Between-study variance: {np.mean(between_vars):.3f} ± {np.std(between_vars):.3f}")
        print(f"   Technical genes: {len(results['technical_genes'])}")
        print(f"   Biological genes: {len(results['biological_genes'])}")
        
        self.results['variance_decomposition'] = results
        return results
    def experiment_2_pathway_activity_profiling(self, dat_lst, label_lst, gene_names):
        """
        Experiment 2: Pathway Activity Profiling
        
        Analyzes pathway activity patterns across studies to understand
        why POSSE's pathway-based trust system is failing
        """
        print("🔬 Experiment 2: Pathway Activity Analysis")
        
        # Define key pathway gene sets
        tb_pathways = {
            'interferon_response': ['IFNG', 'IRF1', 'IRF7', 'STAT1', 'GBP1', 'GBP2', 'OAS1', 'MX1', 'ISG15'],
            'inflammatory_response': ['TNF', 'IL1B', 'IL6', 'NFKB1', 'RELA', 'PTGS2', 'ICAM1', 'CCL2', 'CCL3'],
            'immune_activation': ['CD14', 'CD68', 'TLR2', 'TLR4', 'MYD88', 'IL10', 'IL12A', 'CXCL9', 'CXCL10'],
            'housekeeping': ['ACTB', 'GAPDH', 'B2M', 'HPRT1', 'TBP', 'YWHAZ', 'RPL13A', 'SDHA', 'UBC']
        }
        
        results = {
            'pathway_scores': {},
            'pathway_variance': {},
            'study_signatures': {},
            'tb_vs_control': {}
        }
        
        # Calculate pathway activity scores for each study
        for study, data in dat_lst.items():
            study_results = {}
            
            for pathway_name, pathway_genes in tb_pathways.items():
                # Find genes present in data
                present_genes = [g for g in pathway_genes if g in gene_names]
                
                if len(present_genes) > 0:
                    gene_indices = [np.where(gene_names == g)[0][0] for g in present_genes]
                    pathway_data = data[gene_indices, :]
                    
                    # Calculate pathway activity (mean expression)
                    pathway_activity = np.mean(pathway_data, axis=0)
                    study_results[pathway_name] = pathway_activity
                    
                    print(f"   {study} - {pathway_name}: {len(present_genes)} genes, "
                          f"activity = {np.mean(pathway_activity):.3f} ± {np.std(pathway_activity):.3f}")
            
            results['pathway_scores'][study] = study_results
        
        # Analyze pathway variance across studies
        for pathway_name in tb_pathways.keys():
            pathway_means = []
            pathway_vars = []
            
            for study in dat_lst.keys():
                if pathway_name in results['pathway_scores'][study]:
                    activity = results['pathway_scores'][study][pathway_name]
                    pathway_means.append(np.mean(activity))
                    pathway_vars.append(np.var(activity))
            
            if len(pathway_means) > 1:
                between_study_var = np.var(pathway_means)
                within_study_var = np.mean(pathway_vars)
                
                results['pathway_variance'][pathway_name] = {
                    'between_study': between_study_var,
                    'within_study': within_study_var,
                    'ratio': between_study_var / (within_study_var + 1e-8)
                }
                
                print(f"   {pathway_name} variance ratio (between/within): "
                      f"{results['pathway_variance'][pathway_name]['ratio']:.3f}")
        
        self.results['pathway_analysis'] = results
        return results
    def experiment_3_selective_correction_strategy(self, dat_lst, label_lst, gene_names):
        """
        Experiment 3: Selective Correction Strategy
        
        Tests POSSE variants that preserve biological differences:
        1. Technical-only POSSE (correct only housekeeping/ribosomal genes)
        2. Population-aware POSSE (preserve known population differences)
        3. Hybrid approach (ComBat for technical, POSSE for biological)
        """
        print("🔬 Experiment 3: Selective Correction Strategy")
        
        results = {
            'technical_only': {},
            'population_aware': {},
            'hybrid_approach': {},
            'classification_performance': {}
        }
        
        # Define gene categories
        technical_genes = ['ACTB', 'GAPDH', 'B2M', 'HPRT1', 'TBP', 'YWHAZ', 'RPL13A', 'SDHA', 'UBC']
        ribosomal_genes = [g for g in gene_names if g.startswith('RPS') or g.startswith('RPL')]
        
        # Get study pairs for testing
        study_names = list(dat_lst.keys())
        
        for i, ref_study in enumerate(study_names[:-1]):
            for target_study in study_names[i+1:]:
                print(f"   Testing {ref_study} → {target_study}")
                
                ref_data = BatchData(data=dat_lst[ref_study], gene_indices=gene_names)
                target_data = BatchData(data=dat_lst[target_study], gene_indices=gene_names)
                
                # Strategy 1: Technical-only correction
                tech_genes_present = [g for g in technical_genes + ribosomal_genes if g in gene_names]
                
                if len(tech_genes_present) > 10:
                    # Create restricted pathway dict with only technical genes
                    tech_pathways = {'technical': tech_genes_present}
                    
                    posse_tech = POSSE(
                        pathway_dict=tech_pathways,
                        hyperparams=POSSEHyperparameters(tau=15.0, top_k_percent=0.3)
                    )
                    
                    try:
                        with SuppressOutput():
                            corrected_tech, metadata_tech = posse_tech.align(ref_data, target_data)
                        results['technical_only'][f"{ref_study}_{target_study}"] = {
                            'alpha_mean': metadata_tech.get('alpha_mean', 1.0),
                            'beta_mean': metadata_tech.get('beta_mean', 0.0),
                            'n_genes_corrected': len(tech_genes_present)
                        }
                        print(f"     Technical-only: α={metadata_tech.get('alpha_mean', 1.0):.3f}")
                    except Exception as e:
                        print(f"     Technical-only failed: {e}")
                
                # Strategy 2: Population-aware (preserve immune/inflammatory pathways)
                preserve_pathways = {
                    'interferon': ['IFNG', 'IRF1', 'IRF7', 'STAT1', 'GBP1'],
                    'inflammation': ['TNF', 'IL1B', 'IL6', 'NFKB1', 'RELA']
                }
                
                # Only correct housekeeping genes, preserve immune pathways
                correction_pathways = {'housekeeping': technical_genes}
                
                posse_pop = POSSE(
                    pathway_dict=correction_pathways,
                    hyperparams=POSSEHyperparameters(tau=10.0, top_k_percent=0.4)
                )
                
                try:
                    with SuppressOutput():
                        corrected_pop, metadata_pop = posse_pop.align(ref_data, target_data)
                    results['population_aware'][f"{ref_study}_{target_study}"] = {
                        'alpha_mean': metadata_pop.get('alpha_mean', 1.0),
                        'beta_mean': metadata_pop.get('beta_mean', 0.0)
                    }
                    print(f"     Population-aware: α={metadata_pop.get('alpha_mean', 1.0):.3f}")
                except Exception as e:
                    print(f"     Population-aware failed: {e}")
        
        self.results['selective_correction'] = results
        return results
    def experiment_4_trust_system_debugging(self, dat_lst, label_lst, gene_names):
        """
        Experiment 4: Trust System Debugging
        
        Analyzes POSSE's trust calculations to understand why it's
        over-correcting when housekeeping genes suggest no correction needed
        """
        print("🔬 Experiment 4: Trust System Debugging")
        
        results = {
            'pathway_trust_scores': {},
            'gene_trust_scores': {},
            'trust_vs_biology': {},
            'housekeeping_evidence': {}
        }
        
        # Create a modified POSSE that exposes trust calculations
        class DiagnosticPOSSE(POSSE):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.trust_diagnostics = {}
            
            def pathway_execution(self, X_prime, Y_prime, pathway_indices, C_null, global_gene_idxs, alpha_prior_vec, beta_prior_vec, tau_override=None, top_k_override=None):
                # Call parent method with correct signature
                omega, alpha_est, beta_est, k_raw, metrics = super().pathway_execution(
                    X_prime, Y_prime, pathway_indices, C_null, global_gene_idxs, 
                    alpha_prior_vec, beta_prior_vec, tau_override, top_k_override
                )
                
                # Store trust diagnostics
                pathway_name = f"pathway_{len(self.trust_diagnostics)}"
                self.trust_diagnostics[pathway_name] = {
                    'omega': omega,
                    'alpha_estimates': alpha_est,
                    'beta_estimates': beta_est,
                    'gene_indices': global_gene_idxs,
                    'similarity_scores': k_raw,
                    'metrics': metrics
                }
                
                return omega, alpha_est, beta_est, k_raw, metrics
        
        # Test with diagnostic POSSE
        study_names = list(dat_lst.keys())
        
        if len(study_names) >= 2:
            ref_study = study_names[0]
            target_study = study_names[1]
            
            print(f"   Analyzing trust system: {ref_study} → {target_study}")
            
            ref_data = BatchData(data=dat_lst[ref_study], gene_indices=gene_names)
            target_data = BatchData(data=dat_lst[target_study], gene_indices=gene_names)
            
            # Create pathway dict for analysis
            analysis_pathways = {
                'interferon': ['IFNG', 'IRF1', 'IRF7', 'STAT1', 'GBP1', 'GBP2'],
                'housekeeping': ['ACTB', 'GAPDH', 'B2M', 'HPRT1', 'TBP', 'YWHAZ'],
                'inflammatory': ['TNF', 'IL1B', 'IL6', 'NFKB1', 'RELA', 'PTGS2']
            }
            
            # Filter pathways to genes present in data
            filtered_pathways = {}
            for pathway_name, pathway_genes in analysis_pathways.items():
                present_genes = [g for g in pathway_genes if g in gene_names]
                if len(present_genes) >= 3:
                    filtered_pathways[pathway_name] = present_genes
            
            if len(filtered_pathways) > 0:
                diagnostic_posse = DiagnosticPOSSE(
                    pathway_dict=filtered_pathways,
                    hyperparams=POSSEHyperparameters(tau=25.0, top_k_percent=0.2)
                )
                
                try:
                    with SuppressOutput():
                        corrected_data, metadata = diagnostic_posse.align(ref_data, target_data)
                    
                    # Analyze trust diagnostics
                    for pathway_name, diagnostics in diagnostic_posse.trust_diagnostics.items():
                        results['pathway_trust_scores'][pathway_name] = {
                            'omega': diagnostics['omega'],
                            'mean_alpha': np.mean(diagnostics['alpha_estimates']),
                            'std_alpha': np.std(diagnostics['alpha_estimates']),
                            'mean_similarity': np.mean(diagnostics['similarity_scores'])
                        }
                        
                        print(f"     {pathway_name}: ω={diagnostics['omega']:.3f}, "
                              f"α={np.mean(diagnostics['alpha_estimates']):.3f} ± "
                              f"{np.std(diagnostics['alpha_estimates']):.3f}")
                    
                    # Analyze housekeeping evidence
                    hk_genes = ['ACTB', 'GAPDH', 'B2M', 'HPRT1', 'TBP']
                    hk_present = [g for g in hk_genes if g in gene_names]
                    
                    if len(hk_present) > 0:
                        hk_indices = [np.where(gene_names == g)[0][0] for g in hk_present]
                        
                        ref_hk = ref_data.data[hk_indices, :]
                        target_hk = target_data.data[hk_indices, :]
                        
                        # Calculate what housekeeping genes suggest
                        hk_alpha_suggested = np.std(ref_hk, axis=1) / (np.std(target_hk, axis=1) + 1e-8)
                        hk_beta_suggested = np.mean(ref_hk, axis=1) - hk_alpha_suggested * np.mean(target_hk, axis=1)
                        
                        results['housekeeping_evidence'] = {
                            'suggested_alpha': np.mean(hk_alpha_suggested),
                            'suggested_beta': np.mean(hk_beta_suggested),
                            'alpha_std': np.std(hk_alpha_suggested),
                            'genes_analyzed': hk_present
                        }
                        
                        print(f"     Housekeeping evidence: α={np.mean(hk_alpha_suggested):.3f} ± "
                              f"{np.std(hk_alpha_suggested):.3f}")
                        print(f"     POSSE applied: α={metadata.get('alpha_mean', 1.0):.3f}")
                        
                except Exception as e:
                    print(f"     Trust debugging failed: {e}")
            else:
                print("     No studies available for trust debugging")
        
        self.results['trust_debugging'] = results
        return results
    def experiment_5_population_stratified_validation(self, dat_lst, label_lst, gene_names):
        """
        Experiment 5: Population-Stratified Validation
        
        Tests correction methods within vs across populations to understand
        if POSSE's failure is population-specific
        """
        print("🔬 Experiment 5: Population-Stratified Validation")
        
        results = {
            'within_population': {},
            'cross_population': {},
            'population_signatures': {},
            'correction_effectiveness': {}
        }
        
        # Define population groups (based on study geography/demographics)
        population_groups = {
            'western': ['USA', 'GSE37250_SA'],  # Western populations
            'african': ['Africa', 'GSE37250_M'],  # African populations  
            'asian': ['India', 'GSE39941_M']  # Asian populations
        }
        
        # Filter to available studies
        available_studies = set(dat_lst.keys())
        filtered_groups = {}
        for pop_name, studies in population_groups.items():
            available = [s for s in studies if s in available_studies]
            if len(available) >= 1:
                filtered_groups[pop_name] = available
        
        print(f"   Population groups: {filtered_groups}")
        
        # Test within-population correction
        for pop_name, studies in filtered_groups.items():
            if len(studies) >= 2:
                print(f"   Testing within-population: {pop_name}")
                
                ref_study = studies[0]
                target_study = studies[1]
                
                ref_data = BatchData(data=dat_lst[ref_study], gene_indices=gene_names)
                target_data = BatchData(data=dat_lst[target_study], gene_indices=gene_names)
                
                # Test POSSE within population (suppress output)
                posse_within = POSSE(
                    pathway_dict=self._get_default_pathways(gene_names),
                    hyperparams=POSSEHyperparameters(tau=25.0, top_k_percent=0.2)
                )
                
                try:
                    with SuppressOutput():
                        corrected_within, metadata_within = posse_within.align(ref_data, target_data)
                    
                    results['within_population'][pop_name] = {
                        'alpha_mean': metadata_within.get('alpha_mean', 1.0),
                        'beta_mean': metadata_within.get('beta_mean', 0.0),
                        'ref_study': ref_study,
                        'target_study': target_study
                    }
                    
                    print(f"     Within {pop_name}: α={metadata_within.get('alpha_mean', 1.0):.3f}")
                    
                except Exception as e:
                    print(f"     Within-population correction failed: {e}")
        
        # Test cross-population correction
        pop_names = list(filtered_groups.keys())
        for i, ref_pop in enumerate(pop_names):
            for target_pop in pop_names[i+1:]:
                if len(filtered_groups[ref_pop]) > 0 and len(filtered_groups[target_pop]) > 0:
                    print(f"   Testing cross-population: {ref_pop} → {target_pop}")
                    
                    ref_study = filtered_groups[ref_pop][0]
                    target_study = filtered_groups[target_pop][0]
                    
                    ref_data = BatchData(data=dat_lst[ref_study], gene_indices=gene_names)
                    target_data = BatchData(data=dat_lst[target_study], gene_indices=gene_names)
                    
                    # Test POSSE cross-population (suppress output)
                    posse_cross = POSSE(
                        pathway_dict=self._get_default_pathways(gene_names),
                        hyperparams=POSSEHyperparameters(tau=25.0, top_k_percent=0.2)
                    )
                    
                    try:
                        with SuppressOutput():
                            corrected_cross, metadata_cross = posse_cross.align(ref_data, target_data)
                        
                        results['cross_population'][f"{ref_pop}_{target_pop}"] = {
                            'alpha_mean': metadata_cross.get('alpha_mean', 1.0),
                            'beta_mean': metadata_cross.get('beta_mean', 0.0),
                            'ref_study': ref_study,
                            'target_study': target_study
                        }
                        
                        print(f"     Cross {ref_pop}→{target_pop}: α={metadata_cross.get('alpha_mean', 1.0):.3f}")
                        
                    except Exception as e:
                        print(f"     Cross-population correction failed: {e}")
        
        self.results['population_validation'] = results
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
    def run_all_experiments(self, dat_lst, label_lst, gene_names):
        """Run all diagnostic experiments"""
        print("🚀 Starting POSSE Diagnostic Experiments Suite")
        print("="*60)
        
        # Run experiments
        exp1_results = self.experiment_1_variance_decomposition(dat_lst, label_lst, gene_names)
        print()
        
        exp2_results = self.experiment_2_pathway_activity_profiling(dat_lst, label_lst, gene_names)
        print()
        
        exp3_results = self.experiment_3_selective_correction_strategy(dat_lst, label_lst, gene_names)
        print()
        
        exp4_results = self.experiment_4_trust_system_debugging(dat_lst, label_lst, gene_names)
        print()
        
        exp5_results = self.experiment_5_population_stratified_validation(dat_lst, label_lst, gene_names)
        print()
        
        # Generate summary report
        self.generate_summary_report()
        
        return self.results
    
    def generate_summary_report(self):
        """Generate comprehensive summary of all experiments"""
        print("📊 DIAGNOSTIC SUMMARY REPORT")
        print("="*60)
        
        # Variance decomposition insights
        if 'variance_decomposition' in self.results:
            vd = self.results['variance_decomposition']
            within_vars = list(vd['within_study_var'].values())
            between_vars = list(vd['between_study_var'].values())
            
            print("🔍 VARIANCE DECOMPOSITION:")
            print(f"   Mean within-study variance: {np.mean(within_vars):.3f}")
            print(f"   Mean between-study variance: {np.mean(between_vars):.3f}")
            print(f"   Variance ratio (between/within): {np.mean(between_vars)/np.mean(within_vars):.3f}")
            print(f"   Technical genes identified: {len(vd['technical_genes'])}")
            print(f"   Biological genes identified: {len(vd['biological_genes'])}")
            print()
        
        # Pathway analysis insights
        if 'pathway_analysis' in self.results:
            pa = self.results['pathway_analysis']
            
            print("🧬 PATHWAY ANALYSIS:")
            for pathway, variance_info in pa['pathway_variance'].items():
                ratio = variance_info['ratio']
                print(f"   {pathway}: variance ratio = {ratio:.3f}")
            print()
        
        # Trust system insights
        if 'trust_debugging' in self.results:
            td = self.results['trust_debugging']
            
            print("🎯 TRUST SYSTEM ANALYSIS:")
            if 'housekeeping_evidence' in td and td['housekeeping_evidence']:
                hk = td['housekeeping_evidence']
                if 'suggested_alpha' in hk:
                    print(f"   Housekeeping suggests α = {hk['suggested_alpha']:.3f}")
                    print(f"   Genes analyzed: {hk['genes_analyzed']}")
                else:
                    print("   Housekeeping analysis incomplete")
            
            if 'pathway_trust_scores' in td and td['pathway_trust_scores']:
                print("   Pathway trust scores:")
                for pathway, scores in td['pathway_trust_scores'].items():
                    print(f"     {pathway}: ω={scores['omega']:.3f}, α={scores['mean_alpha']:.3f}")
            else:
                print("   Trust system debugging failed - see experiment 4 output")
            print()
        
        # Population validation insights
        if 'population_validation' in self.results:
            pv = self.results['population_validation']
            
            print("🌍 POPULATION VALIDATION:")
            
            if 'within_population' in pv:
                print("   Within-population corrections:")
                for pop, result in pv['within_population'].items():
                    print(f"     {pop}: α={result['alpha_mean']:.3f}")
            
            if 'cross_population' in pv:
                print("   Cross-population corrections:")
                for comparison, result in pv['cross_population'].items():
                    print(f"     {comparison}: α={result['alpha_mean']:.3f}")
            print()
        
        print("🎯 KEY INSIGHTS FOR POSSE v6.0:")
        print("   1. Check if between-study variance >> within-study variance")
        print("   2. Identify which pathways drive over-correction")
        print("   3. Validate if housekeeping evidence is being ignored")
        print("   4. Test if population differences are being treated as batch effects")
        print("   5. Consider selective correction strategies")

def main():
    """Main function to run diagnostic experiments"""
    
    # Load TB data
    print("Loading TB data...")
    
    # This would load real data - for now using placeholder
    # In practice, this should load from TB_real_data.RData
    dat_lst = {}  # Will be populated with real data
    label_lst = {}  # Will be populated with real labels  
    gene_names = np.array([])  # Will be populated with real gene names
    
    print("Note: This is a framework - needs real TB data to run")
    print("To use: load TB_real_data.RData and populate dat_lst, label_lst, gene_names")
    
    # Initialize diagnostic suite
    config = ExperimentConfig(
        output_dir="posse_diagnostics",
        save_plots=True,
        verbose=True
    )
    
    diagnostic_suite = POSSEDiagnosticSuite(config)
    
    # Run experiments (when real data is available)
    # results = diagnostic_suite.run_all_experiments(dat_lst, label_lst, gene_names)
    
    print("Diagnostic framework ready. Load real data to execute experiments.")

if __name__ == "__main__":
    main()