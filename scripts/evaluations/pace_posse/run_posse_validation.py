#!/usr/bin/env python3
"""
Runner script for POSSE Validation Suite v2.0
Loads real TB data and executes rigorous method validation
"""

import sys
sys.path.append('scripts')
sys.path.append('.')

import numpy as np
import pandas as pd
import subprocess
import pickle
import os
from posse_validation_suite import POSSEValidationSuite, ValidationConfig

def load_tb_data():
    """Load TB data using R script"""
    print("Loading real TB data...")
    
    # Use R to load the data
    r_script = """
    load('data/TB_real_data.RData')
    
    # Print available objects
    cat("Available objects:", ls(), "\\n")
    
    # Print study information
    for(study in names(dat_lst)) {
        cat("Study", study, "- Data shape:", dim(dat_lst[[study]]), "Labels:", length(label_lst[[study]]), "\\n")
    }
    
    # Save gene names
    gene_names <- rownames(dat_lst[[1]])
    write.table(gene_names, 'temp_gene_names.txt', row.names=FALSE, col.names=FALSE, quote=FALSE)
    cat("Saved", length(gene_names), "gene names\\n")
    
    # Save each study's data and labels
    for(study in names(dat_lst)) {
        # Save data matrix
        write.table(dat_lst[[study]], paste0('temp_data_', study, '.txt'), 
                   row.names=FALSE, col.names=FALSE, sep='\\t')
        
        # Save labels
        write.table(label_lst[[study]], paste0('temp_labels_', study, '.txt'), 
                   row.names=FALSE, col.names=FALSE)
        
        cat("Saved data for study:", study, "\\n")
    }
    """
    
    # Execute R script
    result = subprocess.run(['R', '--slave', '--vanilla'], 
                          input=r_script, text=True, capture_output=True)
    
    if result.returncode != 0:
        print(f"R script failed: {result.stderr}")
        return None, None, None
    
    print("R output:", result.stdout)
    
    # Load the saved data back into Python
    dat_lst = {}
    label_lst = {}
    
    # Load gene names
    try:
        gene_names = np.loadtxt('temp_gene_names.txt', dtype=str)
        print(f"Loaded {len(gene_names)} gene names")
    except:
        print("Failed to load gene names")
        return None, None, None
    
    # Find all data files
    data_files = [f for f in os.listdir('.') if f.startswith('temp_data_')]
    
    for data_file in data_files:
        # Extract study name
        study = data_file.replace('temp_data_', '').replace('.txt', '')
        
        try:
            # Load data matrix
            data_matrix = np.loadtxt(data_file, delimiter='\t')
            dat_lst[study] = data_matrix
            
            # Load labels
            label_file = f'temp_labels_{study}.txt'
            if os.path.exists(label_file):
                labels = np.loadtxt(label_file)
                # Handle single value case
                if labels.ndim == 0:
                    labels = np.array([labels])
                label_lst[study] = labels
            else:
                # Create dummy labels if not available
                label_lst[study] = np.zeros(data_matrix.shape[1])
            
            print(f"Loaded {study}: {data_matrix.shape} data, {len(label_lst[study])} labels")
            
        except Exception as e:
            print(f"Failed to load study {study}: {e}")
    
    # Clean up temporary files
    temp_files = [f for f in os.listdir('.') if f.startswith('temp_')]
    for temp_file in temp_files:
        try:
            os.remove(temp_file)
        except:
            pass
    
    return dat_lst, label_lst, gene_names

def main():
    """Main execution function"""
    print("🚀 POSSE Validation Suite v2.0 with Real TB Data")
    print("="*60)
    
    # Load TB data
    dat_lst, label_lst, gene_names = load_tb_data()
    
    if dat_lst is None:
        print("❌ Failed to load TB data")
        return
    
    print(f"✅ Loaded data for {len(dat_lst)} studies")
    print(f"✅ Gene names: {len(gene_names)}")
    
    # Initialize validation suite with strict thresholds
    config = ValidationConfig(
        output_dir="posse_validation_results",
        save_plots=True,
        verbose=True,
        # Strict validation thresholds
        trust_threshold=0.5,
        alpha_deviation_threshold=0.2,
        silence_rmse_threshold=1e-9,
        topology_preservation_threshold=0.1
    )
    
    validation_suite = POSSEValidationSuite(config)
    
    # Run complete validation suite
    results = validation_suite.run_validation_suite(dat_lst, label_lst, gene_names)
    
    # Save results
    os.makedirs(config.output_dir, exist_ok=True)
    
    # Save detailed results
    with open(f"{config.output_dir}/validation_results.pkl", 'wb') as f:
        pickle.dump(results, f)
    
    # Save failure summary
    with open(f"{config.output_dir}/validation_failures.txt", 'w') as f:
        f.write("POSSE Validation Suite v2.0 - Failure Report\n")
        f.write("="*50 + "\n\n")
        
        if len(validation_suite.failures) == 0:
            f.write("🎉 ALL TESTS PASSED - No validation failures detected!\n")
        else:
            f.write(f"❌ {len(validation_suite.failures)} VALIDATION FAILURES:\n\n")
            for i, failure in enumerate(validation_suite.failures, 1):
                f.write(f"{i}. {failure}\n")
    
    print(f"✅ Validation completed!")
    print(f"📁 Results saved to: {config.output_dir}")
    print(f"💾 Detailed results saved to: {config.output_dir}/validation_results.pkl")
    print(f"📄 Failure report saved to: {config.output_dir}/validation_failures.txt")

if __name__ == "__main__":
    main()