#!/usr/bin/env python3
"""
Run POSSE Diagnostic Experiments with Real TB Data
"""

import sys
sys.path.append('scripts')

import subprocess
import tempfile
import os
import numpy as np
import pandas as pd
from posse_diagnostic_experiments import POSSEDiagnosticSuite, ExperimentConfig

def load_real_tb_data():
    """Load real TB data from the TB_real_data.RData file"""
    
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
    
    return dat_lst, label_lst, gene_names

def main():
    """Main function to run POSSE diagnostics with real data"""
    
    print("🚀 POSSE Diagnostic Experiments with Real TB Data")
    print("="*60)
    
    # Load real TB data
    try:
        dat_lst, label_lst, gene_names = load_real_tb_data()
        
        if len(dat_lst) == 0:
            print("❌ No data loaded. Check if TB_real_data.RData exists in data/ directory")
            return
        
        print(f"✅ Loaded data for {len(dat_lst)} studies")
        print(f"✅ Gene names: {len(gene_names) if gene_names is not None else 0}")
        
    except Exception as e:
        print(f"❌ Failed to load TB data: {e}")
        return
    
    # Initialize diagnostic suite
    config = ExperimentConfig(
        output_dir="posse_diagnostics_results",
        save_plots=True,
        verbose=True
    )
    
    diagnostic_suite = POSSEDiagnosticSuite(config)
    
    # Run all experiments
    try:
        results = diagnostic_suite.run_all_experiments(dat_lst, label_lst, gene_names)
        
        print("✅ All experiments completed successfully!")
        print(f"📁 Results saved to: {config.output_dir}")
        
        # Save results to file
        import pickle
        results_file = os.path.join(config.output_dir, "diagnostic_results.pkl")
        os.makedirs(config.output_dir, exist_ok=True)
        
        with open(results_file, 'wb') as f:
            pickle.dump(results, f)
        
        print(f"💾 Detailed results saved to: {results_file}")
        
    except Exception as e:
        print(f"❌ Experiments failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()