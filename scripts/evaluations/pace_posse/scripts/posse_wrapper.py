#!/usr/bin/env python3
"""
POSSE wrapper for R integration
Provides a simple interface to call POSSE from R via reticulate
Updated to use POSSE v6.0 with gene-wise standardization
"""

import numpy as np
import sys
import os

# Add the scripts directory to path to import posse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from posse_v6 import POSSEv6, BatchData, POSSEv6Hyperparameters

def check_and_normalize(ref_data, target_data):
    """
    Safety Valve: Detects and fixes massive library size mismatches (e.g., 381x).
    POSSE requires data to be roughly comparable (linear scale) to find peers.
    """
    # 1. Calculate Global Scale Difference
    # Use median of non-zeros to avoid sparse noise
    ref_med = np.median(ref_data.data[ref_data.data > 0])
    tgt_med = np.median(target_data.data[target_data.data > 0])
    
    if tgt_med == 0 or ref_med == 0:
        return ref_data, target_data # Cannot compute, abort safety check
        
    ratio = ref_med / tgt_med
    
    print(f"  [Safety Check] Global Scale Ratio: {ratio:.2f}x")
    
    # 2. Trigger Normalization if > 10x difference
    # This implies one is Raw Counts and one is Normalized, or massive depth diff.
    if ratio > 10.0 or ratio < 0.1:
        print("  [WARNING] Massive scale difference detected! Auto-converting to CPM.")
        
        def to_cpm(data):
            # Sum per sample
            lib_size = np.sum(data, axis=0)
            # Avoid div/0
            lib_size[lib_size == 0] = 1.0
            # Convert to CPM (scale to 1 million)
            return data / lib_size[np.newaxis, :] * 1e6
            
        # Apply to both to ensure same scale
        # We modify the .data attribute in place or create new objects
        ref_data.data = to_cpm(ref_data.data)
        target_data.data = to_cpm(target_data.data)
        
        # Re-check
        new_ref = np.median(ref_data.data[ref_data.data > 0])
        new_tgt = np.median(target_data.data[target_data.data > 0])
        print(f"  [FIXED] New Scale Ratio: {new_ref/new_tgt:.2f}x")
        
    return ref_data, target_data

def get_posse_config(posse_variant):
    """
    Get POSSE v6.0 configuration parameters based on variant name.
    v6.0 uses gene-wise standardization for robust cross-study correction.
    
    Args:
        posse_variant: str - one of 'default', 'aggressive', 'focused', 'conservative', 'housekeeping'
    
    Returns:
        dict with POSSE configuration parameters
    """
    configs = {
        'default': {
            'pathway_source': 'MSigDB_Hallmark_2020',
            'tau': 25.0,
            'top_k_percent': 0.20,
            'max_iter': 3,
            'min_pathway_size': 5
        },
        'aggressive': {
            'pathway_source': 'MSigDB_Hallmark_2020',
            'tau': 50.0,
            'top_k_percent': 0.10,
            'max_iter': 4,
            'min_pathway_size': 5
        },
        'focused': {
            'pathway_source': 'MSigDB_Hallmark_2020',
            'tau': 30.0,
            'top_k_percent': 0.15,
            'max_iter': 3,
            'min_pathway_size': 5
        },
        'conservative': {
            'pathway_source': 'MSigDB_Hallmark_2020',
            'tau': 15.0,
            'top_k_percent': 0.25,
            'max_iter': 2,
            'min_pathway_size': 5
        },
        'housekeeping': {
            'pathway_source': 'MSigDB_Hallmark_2020',
            'tau': 20.0,
            'top_k_percent': 0.20,
            'max_iter': 3,
            'min_pathway_size': 5
        },
        'ultra_aggressive': {
            'pathway_source': 'MSigDB_Hallmark_2020',
            'tau': 75.0,
            'top_k_percent': 0.05,
            'max_iter': 5,
            'min_pathway_size': 5
        }
    }
    
    return configs.get(posse_variant, configs['default'])

def apply_posse_correction(train_data, test_data, train_batch, posse_variant='default', gene_names=None, save_diagnostics=None):
    """
    Apply POSSE v6.0 batch correction to training and test data
    
    Args:
        train_data: numpy array (genes x samples) - training data
        test_data: numpy array (genes x samples) - test data  
        train_batch: numpy array - batch labels for training samples
        posse_variant: str - POSSE variant ('default', 'aggressive', 'focused', 'conservative', 'housekeeping', 'ultra_aggressive')
        gene_names: list/array - gene names corresponding to rows (optional)
        save_diagnostics: str - path to save diagnostic CSV file (optional)
    
    Returns:
        dict with 'train_corrected' and 'test_corrected' arrays
    """
    
    # Convert to numpy arrays if needed
    train_data = np.asarray(train_data, dtype=np.float64)
    test_data = np.asarray(test_data, dtype=np.float64)
    train_batch = np.asarray(train_batch, dtype=np.int32)
    
    print(f"POSSE v6.0 {posse_variant} input shapes: train {train_data.shape}, test {test_data.shape}")
    print(f"Train batch unique values: {np.unique(train_batch)}")
    
    # Get configuration for this POSSE variant
    config = get_posse_config(posse_variant)
    
    # Create gene indices - use provided gene names or row indices as strings
    n_genes = train_data.shape[0]
    if gene_names is not None:
        gene_indices = np.asarray(gene_names, dtype=str)
        print(f"Using provided gene names. Sample genes: {gene_indices[:3]}")
    else:
        gene_indices = np.arange(n_genes).astype(str)
        print(f"No gene names provided, using row indices: {gene_indices[:3]}")
    
    if len(gene_indices) != n_genes:
        print(f"Warning: Gene names length ({len(gene_indices)}) doesn't match data rows ({n_genes})")
        gene_indices = np.arange(n_genes).astype(str)
    
    try:
        # Load pathway dictionary
        pathway_dict = load_pathway_dict(config['pathway_source'], gene_indices)
        
        # Initialize POSSE v6.0 with variant-specific parameters
        hyperparams = POSSEv6Hyperparameters(
            tau=config['tau'],
            top_k_percent=config['top_k_percent'],
            max_iter=config['max_iter'],
            min_pathway_size=config['min_pathway_size']
        )
        
        posse_model = POSSEv6(
            pathway_dict=pathway_dict,
            hyperparams=hyperparams
        )
        
        print(f"POSSE v6.0 {posse_variant} initialized with {len(pathway_dict)} pathways")
        print(f"Config: tau={config['tau']}, top_k_percent={config['top_k_percent']}")
        print(f"POSSE v6.0: Gene-wise standardization for robust cross-study correction")
        
        # Use first batch as reference, combine others with test as target
        unique_batches = np.unique(train_batch)
        ref_batch_id = unique_batches[0]
        
        # Reference data (first batch)
        ref_mask = train_batch == ref_batch_id
        ref_data = BatchData(
            data=train_data[:, ref_mask],
            gene_indices=gene_indices
        )
        
        # Target data (remaining batches + test data)
        target_mask = train_batch != ref_batch_id
        if np.any(target_mask):
            target_train = train_data[:, target_mask]
            target_data_matrix = np.hstack([target_train, test_data])
            n_target_train = target_train.shape[1]
        else:
            target_data_matrix = test_data
            n_target_train = 0
            
        target_data = BatchData(
            data=target_data_matrix,
            gene_indices=gene_indices
        )
        
        # Safety check for massive scale differences
        print("Running Pre-flight Checks...")
        ref_data, target_data = check_and_normalize(ref_data, target_data)
        
        print(f"Reference batch samples: {ref_data.data.shape[1]}")
        print(f"Target data samples: {target_data.data.shape[1]} (train: {n_target_train}, test: {test_data.shape[1]})")
        
        # Apply POSSE v6.0 correction
        corrected_batch, metadata = posse_model.align(ref_data, target_data)
        corrected_data = corrected_batch.data
        
        print(f"POSSE v6.0 {posse_variant} correction completed. Output shape: {corrected_data.shape}")
        print(f"Alpha mean: {metadata.get('alpha_mean', 'N/A')}")
        print(f"Beta mean: {metadata.get('beta_mean', 'N/A')}")
        print(f"Genes covered: {metadata.get('genes_covered', 'N/A')}/{metadata.get('genes_total', 'N/A')}")
        
        # Save diagnostics if requested
        if save_diagnostics:
            try:
                diagnostic_data = [{
                    'gene_id': 'summary',
                    'gene_index': -1,
                    'alpha_final': metadata.get('alpha_mean', 1.0),
                    'beta_final': metadata.get('beta_mean', 0.0),
                    'alpha_internal': metadata.get('alpha_internal', 1.0),
                    'beta_internal': metadata.get('beta_internal', 0.0),
                    'variant': posse_variant,
                    'tau': config['tau'],
                    'top_k_percent': config['top_k_percent'],
                    'genes_covered': metadata.get('genes_covered', 0),
                    'genes_total': metadata.get('genes_total', 0),
                    'version': metadata.get('version', '6.0')
                }]
                
                # Write CSV
                import csv
                os.makedirs(os.path.dirname(save_diagnostics), exist_ok=True)
                
                with open(save_diagnostics, 'w', newline='') as f:
                    writer = csv.DictWriter(f, fieldnames=diagnostic_data[0].keys())
                    writer.writeheader()
                    writer.writerows(diagnostic_data)
                
                print(f"✓ Diagnostic file created: {save_diagnostics}")
                
            except Exception as e:
                print(f"WARNING: Failed to save diagnostics: {str(e)}")
        
        # Add variant info to metadata
        metadata['posse_variant'] = posse_variant
        metadata['config'] = config
        
        # Reconstruct training and test data
        train_corrected = np.zeros_like(train_data)
        train_corrected[:, ref_mask] = ref_data.data  # Reference batch unchanged
        
        if n_target_train > 0:
            train_corrected[:, target_mask] = corrected_data[:, :n_target_train]
            test_corrected = corrected_data[:, n_target_train:]
        else:
            test_corrected = corrected_data
        
        return {
            'train_corrected': train_corrected,
            'test_corrected': test_corrected,
            'metadata': metadata
        }
        
    except Exception as e:
        print(f"POSSE v6.0 {posse_variant} correction failed: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        # Return original data if POSSE fails
        return {
            'train_corrected': train_data,
            'test_corrected': test_data,
            'metadata': {'error': str(e), 'posse_variant': posse_variant}
        }


def load_pathway_dict(pathway_source, gene_names):
    """Load pathway dictionary from source, filtered for available genes"""
    # Default pathways if source not available
    default_pathways = {
        'interferon': ['IFNG', 'IRF1', 'IRF7', 'STAT1', 'GBP1', 'GBP2', 'OAS1', 'MX1', 'ISG15', 'IFIT1'],
        'inflammatory': ['TNF', 'IL1B', 'IL6', 'NFKB1', 'RELA', 'PTGS2', 'ICAM1', 'CCL2', 'CXCL8'],
        'housekeeping': ['ACTB', 'GAPDH', 'B2M', 'HPRT1', 'TBP', 'YWHAZ', 'RPL13A', 'SDHA'],
        'apoptosis': ['BCL2', 'BAX', 'CASP3', 'CASP8', 'CASP9', 'TP53', 'FAS', 'FASLG'],
        'cell_cycle': ['CCND1', 'CCNE1', 'CDK2', 'CDK4', 'RB1', 'E2F1', 'MYC', 'PCNA'],
        'metabolism': ['HK1', 'HK2', 'PKM', 'LDHA', 'PFKFB3', 'SLC2A1', 'PDK1'],
        'immune_response': ['CD4', 'CD8A', 'CD19', 'CD14', 'FCGR1A', 'TLR2', 'TLR4', 'IL10']
    }
    
    # Try to load MSigDB pathways if available
    try:
        msigdb_path = os.path.join(os.path.dirname(__file__), 'msigdb_hallmark.json')
        if os.path.exists(msigdb_path):
            import json
            with open(msigdb_path, 'r') as f:
                pathways = json.load(f)
            print(f"Loaded {len(pathways)} pathways from MSigDB")
        else:
            pathways = default_pathways
            print(f"Using {len(pathways)} default pathways")
    except:
        pathways = default_pathways
        print(f"Using {len(pathways)} default pathways")
    
    # Filter pathways to only include genes present in the data
    gene_set = set(gene_names)
    filtered_pathways = {}
    for name, genes in pathways.items():
        present_genes = [g for g in genes if g in gene_set]
        if len(present_genes) >= 3:
            filtered_pathways[name] = present_genes
    
    print(f"Filtered to {len(filtered_pathways)} pathways with ≥3 genes present")
    return filtered_pathways

if __name__ == "__main__":
    # Test the wrapper
    print("POSSE wrapper loaded successfully")