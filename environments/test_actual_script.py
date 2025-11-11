#!/usr/bin/env python3
"""
Test actual pipeline functionality with the uv Python environment
"""

import sys
import os
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd

# Add scripts to path
sys.path.insert(0, 'scripts')

def test_data_processing():
    """Test data processing similar to pipeline scripts"""
    print("=" * 60)
    print("Testing Data Processing")
    print("=" * 60)
    
    # Create synthetic gene expression data
    n_samples = 100
    n_genes = 50
    
    df = pd.DataFrame({
        'meta_er_status': np.random.choice(['Positive', 'Negative'], n_samples),
        'meta_batch': np.random.choice(['A', 'B', 'C'], n_samples),
        'meta_source': np.random.choice(['Dataset1', 'Dataset2'], n_samples),
    })
    
    # Add gene expression columns
    for i in range(n_genes):
        df[f'gene_{i}'] = np.random.randn(n_samples)
    
    print(f"✓ Created synthetic dataset: {df.shape}")
    
    # Test filtering
    positive_samples = df[df['meta_er_status'] == 'Positive']
    print(f"✓ Filtered positive samples: {positive_samples.shape[0]}")
    
    # Test grouping (only numeric columns)
    gene_cols = [col for col in df.columns if col.startswith('gene_')]
    grouped = df.groupby('meta_batch')[gene_cols].mean()
    print(f"✓ Grouped by batch: {grouped.shape}")
    
    # Test CSV I/O
    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / 'test_data.csv'
        df.to_csv(test_file, index=False)
        df_loaded = pd.read_csv(test_file)
        assert df_loaded.shape == df.shape
        print(f"✓ CSV I/O successful")
    
    return True

def test_machine_learning():
    """Test machine learning similar to pipeline scripts"""
    print("\n" + "=" * 60)
    print("Testing Machine Learning")
    print("=" * 60)
    
    from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
    from sklearn.model_selection import cross_val_score
    from sklearn.metrics import accuracy_score, roc_auc_score
    
    # Create synthetic data
    n_samples = 200
    n_features = 50
    
    X = np.random.randn(n_samples, n_features)
    y = np.random.randint(0, 2, n_samples)
    
    print(f"✓ Created synthetic ML dataset: X={X.shape}, y={y.shape}")
    
    # Test Random Forest (used in pipeline)
    rf = RandomForestClassifier(n_estimators=50, random_state=42)
    rf_scores = cross_val_score(rf, X, y, cv=3, scoring='accuracy')
    print(f"✓ Random Forest CV accuracy: {rf_scores.mean():.3f} ± {rf_scores.std():.3f}")
    
    # Test Histogram Gradient Boosting (used in pipeline)
    hgb = HistGradientBoostingClassifier(max_iter=50, random_state=42)
    hgb_scores = cross_val_score(hgb, X, y, cv=3, scoring='accuracy')
    print(f"✓ HistGradientBoosting CV accuracy: {hgb_scores.mean():.3f} ± {hgb_scores.std():.3f}")
    
    # Test ROC AUC (used in pipeline)
    rf.fit(X, y)
    y_pred_proba = rf.predict_proba(X)[:, 1]
    roc_auc = roc_auc_score(y, y_pred_proba)
    print(f"✓ ROC AUC score: {roc_auc:.3f}")
    
    return True

def test_local_modules():
    """Test importing local pipeline modules"""
    print("\n" + "=" * 60)
    print("Testing Local Module Imports")
    print("=" * 60)
    
    try:
        from utils import DataFrameCache, HashCache
        print("✓ Imported utils.DataFrameCache")
        print("✓ Imported utils.HashCache")
        
        # Test instantiation
        cache = DataFrameCache()
        print("✓ Created DataFrameCache instance")
        
        from evaluations.util import repeated_cross_val
        print("✓ Imported evaluations.util.repeated_cross_val")
        
        return True
    except ImportError as e:
        print(f"✗ Failed to import local modules: {e}")
        return False

def test_dimensionality_reduction():
    """Test dimensionality reduction (used in pipeline)"""
    print("\n" + "=" * 60)
    print("Testing Dimensionality Reduction")
    print("=" * 60)
    
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    import umap
    
    # Create high-dimensional data
    X = np.random.randn(100, 50)
    
    # Test PCA
    pca = PCA(n_components=10)
    X_pca = pca.fit_transform(X)
    print(f"✓ PCA: {X.shape} → {X_pca.shape}")
    
    # Test t-SNE
    tsne = TSNE(n_components=2, random_state=42)
    X_tsne = tsne.fit_transform(X)
    print(f"✓ t-SNE: {X.shape} → {X_tsne.shape}")
    
    # Test UMAP
    reducer = umap.UMAP(n_components=2, random_state=42)
    X_umap = reducer.fit_transform(X)
    print(f"✓ UMAP: {X.shape} → {X_umap.shape}")
    
    return True

def test_pytorch():
    """Test PyTorch functionality (used in pipeline)"""
    print("\n" + "=" * 60)
    print("Testing PyTorch")
    print("=" * 60)
    
    import torch
    import torch.nn as nn
    
    # Create simple neural network
    model = nn.Sequential(
        nn.Linear(50, 128),
        nn.ReLU(),
        nn.Dropout(0.2),
        nn.Linear(128, 64),
        nn.ReLU(),
        nn.Linear(64, 2)
    )
    
    print(f"✓ Created neural network with {sum(p.numel() for p in model.parameters())} parameters")
    
    # Test forward pass
    x = torch.randn(32, 50)
    y = model(x)
    print(f"✓ Forward pass: {x.shape} → {y.shape}")
    
    # Check CUDA
    cuda_available = torch.cuda.is_available()
    print(f"✓ CUDA available: {cuda_available}")
    
    return True

def main():
    """Run all tests"""
    print("\n" + "=" * 60)
    print("Pipeline Functionality Test")
    print("=" * 60)
    print(f"Python: {sys.version}")
    print(f"Working directory: {os.getcwd()}")
    print()
    
    tests = [
        ("Data Processing", test_data_processing),
        ("Machine Learning", test_machine_learning),
        ("Local Modules", test_local_modules),
        ("Dimensionality Reduction", test_dimensionality_reduction),
        ("PyTorch", test_pytorch),
    ]
    
    all_passed = True
    for test_name, test_func in tests:
        try:
            if not test_func():
                all_passed = False
                print(f"\n✗ {test_name} test FAILED")
        except Exception as e:
            all_passed = False
            print(f"\n✗ {test_name} test FAILED with exception:")
            print(f"  {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
    
    # Summary
    print("\n" + "=" * 60)
    if all_passed:
        print("✓ ALL PIPELINE FUNCTIONALITY TESTS PASSED")
        print("=" * 60)
        print("\nThe Python environment is fully functional for pipeline use!")
        return 0
    else:
        print("✗ SOME TESTS FAILED")
        print("=" * 60)
        return 1

if __name__ == "__main__":
    sys.exit(main())
