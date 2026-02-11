#!/usr/bin/env python3
"""
PACE wrapper for R integration
Provides a simple interface to call PACE from R via reticulate
"""

import numpy as np
import sys
import os

# Add the scripts directory to path to import pace
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pace import PACE_v22, BatchData, PACEHyperparameters

def get_pace_config(pace_variant):
    """
    Get PACE configuration parameters based on variant name
    
    Args:
        pace_variant: str - one of 'default', 'aggressive', 'focused', 'conservative'
    
    Returns:
        dict with PACE configuration parameters
    """
    configs = {
        'default': {
            'pathway_source': 'data/pathways/comprehensive.gmt',
            'w_prior': 1.0,
            'tau': 10.0,
            'use_topk_masking': False,  # Keep soft attention for default
            'topk_ratio': 0.20
        },
        'aggressive': {
            'pathway_source': 'data/pathways/comprehensive.gmt',
            'w_prior': 2.0,  # CORRECTED: Higher prior for safety net with high tau
            'tau': 50.0,     # High contrast for aggressive matching
            'use_topk_masking': True,   # Enable Top-K masking to fix Simpson's Paradox
            'topk_ratio': 0.10,  # Keep only top 10% of reference samples (more aggressive)
            'similarity_metric': 'centered_cosine'  # Use Pearson correlation for better biological matching
        },
        'focused': {
            'pathway_source': 'data/pathways/comprehensive.gmt',  # Use comprehensive pathways for better gene coverage
            'w_prior': 1.5,  # CORRECTED: Moderate prior for moderate tau
            'tau': 30.0,     # Medium-high contrast for focused matching
            'use_topk_masking': True,   # Enable Top-K masking
            'topk_ratio': 0.15,  # Keep top 15% of reference samples
            'similarity_metric': 'centered_cosine'  # Use Pearson correlation
        },
        'conservative': {
            'pathway_source': 'data/pathways/comprehensive.gmt',
            'w_prior': 0.5,  # Lower prior for conservative approach (less regularization needed)
            'tau': 5.0,      # Lower contrast for conservative approach
            'use_topk_masking': False,  # Keep soft attention for conservative
            'topk_ratio': 0.30,
            'similarity_metric': 'centered_cosine'
        },
        'aggressive_cosine': {
            'pathway_source': 'data/pathways/comprehensive.gmt',
            'w_prior': 2.0,
            'tau': 50.0,
            'use_topk_masking': True,
            'topk_ratio': 0.10,
            'similarity_metric': 'cosine'  # Test standard cosine similarity
        },
        'focused_cosine': {
            'pathway_source': 'data/pathways/comprehensive.gmt',
            'w_prior': 1.5,
            'tau': 30.0,
            'use_topk_masking': True,
            'topk_ratio': 0.15,
            'similarity_metric': 'cosine'  # Test standard cosine similarity
        },
        'ultra_aggressive': {
            'pathway_source': 'data/pathways/hallmark.gmt',  # Use focused hallmark pathways
            'w_prior': 0.0  # Completely disable global blending
        },
        'extreme_aggressive': {
            'pathway_source': 'data/pathways/hallmark.gmt',  # Focused TB-relevant pathways
            'w_prior': 0.0     # Complete trust in local data
        },
        'iterative_aggressive': {
            'pathway_source': 'data/pathways/hallmark.gmt',  # Focused pathways
            'w_prior': 0.0    # Pure local variance
        }
    }
    
    return configs.get(pace_variant, configs['default'])

def apply_pace_correction(train_data, test_data, train_batch, pace_variant='default', gene_names=None, save_diagnostics=None):
    """
    Apply PACE batch correction to training and test data
    
    Args:
        train_data: numpy array (genes x samples) - training data
        test_data: numpy array (genes x samples) - test data  
        train_batch: numpy array - batch labels for training samples
        pace_variant: str - PACE variant ('default', 'aggressive', 'focused', 'conservative')
        gene_names: list/array - gene names corresponding to rows (optional)
        save_diagnostics: str - path to save diagnostic CSV file (optional)
    
    Returns:
        dict with 'train_corrected' and 'test_corrected' arrays
    """
    
    # Convert to numpy arrays if needed
    train_data = np.asarray(train_data, dtype=np.float64)
    test_data = np.asarray(test_data, dtype=np.float64)
    train_batch = np.asarray(train_batch, dtype=np.int32)
    
    print(f"PACE {pace_variant} input shapes: train {train_data.shape}, test {test_data.shape}")
    print(f"Train batch unique values: {np.unique(train_batch)}")
    
    # Get configuration for this PACE variant
    config = get_pace_config(pace_variant)
    
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
        # Initialize PACE with variant-specific parameters
        hyperparams = PACEHyperparameters(
            w_prior=config['w_prior'],
            tau=config['tau'],
            use_hard_gating=config.get('use_topk_masking', True),
            hard_gating_ratio=config.get('topk_ratio', 0.20),
            similarity_metric=config.get('similarity_metric', 'centered_cosine')
        )
        
        pace_kwargs = {
            'pathway_source': config['pathway_source'],
            'organism': 'Human',
            'hyperparams': hyperparams
        }
        
        pace = PACE_v22(**pace_kwargs)
        
        print(f"PACE {pace_variant} initialized with {len(pace.pathway_dict)} pathways")
        print(f"Config: w_prior={config['w_prior']}, tau={config['tau']}")
        
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
        
        print(f"Reference batch samples: {ref_data.data.shape[1]}")
        print(f"Target data samples: {target_data.data.shape[1]} (train: {n_target_train}, test: {test_data.shape[1]})")
        
        # Apply PACE correction with diagnostic output
        corrected_data, metrics = pace.align(ref_data, target_data, diagnostic_output_path=save_diagnostics)
        
        # Validate that diagnostic file was created if requested - FAIL FAST
        if save_diagnostics:
            if not os.path.exists(save_diagnostics):
                raise RuntimeError(f"CRITICAL: PACE diagnostic file was not created at expected path: {save_diagnostics}")
            
            # Check file is not empty
            file_size = os.path.getsize(save_diagnostics)
            if file_size == 0:
                raise RuntimeError(f"CRITICAL: PACE diagnostic file is empty: {save_diagnostics}")
            
            # Verify file has valid CSV structure
            try:
                with open(save_diagnostics, 'r') as f:
                    first_line = f.readline().strip()
                    if not first_line or ',' not in first_line:
                        raise RuntimeError(f"CRITICAL: PACE diagnostic file has invalid CSV format: {save_diagnostics}")
            except Exception as e:
                raise RuntimeError(f"CRITICAL: Cannot read PACE diagnostic file {save_diagnostics}: {str(e)}")
            
            print(f"✓ Diagnostic file successfully created and validated: {save_diagnostics} ({file_size} bytes)")
        else:
            print("WARNING: No diagnostic output path provided - diagnostics will not be saved")
        
        print(f"PACE {pace_variant} correction completed. Output shape: {corrected_data.shape}")
        print(f"Iterations: {metrics.get('iters', 'N/A')}, Prior strength: {metrics.get('prior_strength', 'N/A')}")
        
        # Add variant info to metrics
        metrics['pace_variant'] = pace_variant
        metrics['config'] = config
        
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
            'metrics': metrics
        }
        
    except Exception as e:
        print(f"PACE {pace_variant} correction failed: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        # Return original data if PACE fails
        return {
            'train_corrected': train_data,
            'test_corrected': test_data,
            'metrics': {'error': str(e), 'pace_variant': pace_variant}
        }

if __name__ == "__main__":
    # Test the wrapper
    print("PACE wrapper loaded successfully")