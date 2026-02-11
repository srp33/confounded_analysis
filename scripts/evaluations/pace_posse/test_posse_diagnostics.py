#!/usr/bin/env python3
"""
Test POSSE Diagnostic Framework with Synthetic Data
"""

import sys
sys.path.append('scripts')

import numpy as np
from posse_diagnostic_experiments import POSSEDiagnosticSuite, ExperimentConfig

def generate_synthetic_tb_data():
    """Generate synthetic TB data for testing the diagnostic framework"""
    
    np.random.seed(42)
    
    # Define synthetic studies
    studies = ['USA', 'Africa', 'India']
    n_genes = 500
    n_samples_per_study = 60
    
    # Generate gene names (mix of real and synthetic)
    real_genes = [
        'IFNG', 'IRF1', 'IRF7', 'STAT1', 'GBP1', 'GBP2', 'OAS1', 'MX1', 'ISG15',
        'TNF', 'IL1B', 'IL6', 'NFKB1', 'RELA', 'PTGS2', 'ICAM1', 'CCL2', 'CCL3',
        'CD14', 'CD68', 'TLR2', 'TLR4', 'MYD88', 'IL10', 'IL12A', 'CXCL9', 'CXCL10',
        'ACTB', 'GAPDH', 'B2M', 'HPRT1', 'TBP', 'YWHAZ', 'RPL13A', 'SDHA', 'UBC'
    ]
    
    synthetic_genes = [f"GENE_{i:04d}" for i in range(n_genes - len(real_genes))]
    gene_names = np.array(real_genes + synthetic_genes)
    
    dat_lst = {}
    label_lst = {}
    
    # Generate data for each study
    for i, study in enumerate(studies):
        # Create TB vs Control labels (40% TB, 60% Control)
        n_tb = int(n_samples_per_study * 0.4)
        n_control = n_samples_per_study - n_tb
        labels = np.array([1] * n_tb + [0] * n_control)
        np.random.shuffle(labels)
        
        # Generate expression data
        baseline_expr = np.random.uniform(5, 12, (n_genes, n_samples_per_study))
        
        # Add TB signature to TB samples
        tb_mask = labels == 1
        tb_genes_idx = [j for j, gene in enumerate(gene_names) if gene in [
            'IFNG', 'IRF1', 'TNF', 'IL1B', 'IL6', 'GBP1', 'STAT1'
        ]]
        
        for gene_idx in tb_genes_idx:
            baseline_expr[gene_idx, tb_mask] += np.random.normal(2.0, 0.5, np.sum(tb_mask))
        
        # Add study-specific batch effects
        if study == 'USA':
            baseline_expr += np.random.normal(0.2, 0.1, baseline_expr.shape)
        elif study == 'Africa':
            baseline_expr += np.random.normal(-0.1, 0.15, baseline_expr.shape)
        elif study == 'India':
            baseline_expr += np.random.normal(0.0, 0.2, baseline_expr.shape)
        
        # Add population-specific signatures
        pop_genes_idx = [j for j, gene in enumerate(gene_names) if gene.startswith('GENE_')][:20]
        
        if study == 'Africa':
            # African population signature
            for gene_idx in pop_genes_idx[:10]:
                baseline_expr[gene_idx, :] += np.random.normal(1.0, 0.3, n_samples_per_study)
        elif study == 'India':
            # Asian population signature  
            for gene_idx in pop_genes_idx[10:]:
                baseline_expr[gene_idx, :] += np.random.normal(0.8, 0.2, n_samples_per_study)
        
        # Convert to linear scale (simulate log-transformed data)
        dat_lst[study] = 2**baseline_expr
        label_lst[study] = labels
    
    return dat_lst, label_lst, gene_names

def main():
    """Test the diagnostic framework with synthetic data"""
    
    print("🧪 Testing POSSE Diagnostic Framework with Synthetic Data")
    print("="*60)
    
    # Generate synthetic data
    dat_lst, label_lst, gene_names = generate_synthetic_tb_data()
    
    print(f"Generated synthetic data:")
    for study, data in dat_lst.items():
        print(f"  {study}: {data.shape} (genes x samples)")
        print(f"    TB cases: {np.sum(label_lst[study])}/{len(label_lst[study])}")
    
    print(f"  Total genes: {len(gene_names)}")
    print()
    
    # Initialize diagnostic suite
    config = ExperimentConfig(
        output_dir="test_posse_diagnostics",
        save_plots=False,  # Don't save plots for test
        verbose=True
    )
    
    diagnostic_suite = POSSEDiagnosticSuite(config)
    
    # Run experiments
    try:
        print("Running diagnostic experiments...")
        results = diagnostic_suite.run_all_experiments(dat_lst, label_lst, gene_names)
        
        print("✅ Test completed successfully!")
        print("🎯 Framework is ready for real TB data")
        
    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()