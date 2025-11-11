#!/usr/bin/env python3
"""
Python Package Import Validation Script
Tests all required packages for the batch correction pipeline
Designed to run on both login and compute nodes

Usage:
    python environments/validate_python_imports.py
    ./environments/run_with_env.sh --sbatch environments/validate_python_imports.py
"""

import sys
import os
import platform
from datetime import datetime
from pathlib import Path

def print_header(title):
    """Print a formatted header"""
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70)

def print_section(title):
    """Print a section header"""
    print(f"\n--- {title} ---")

def test_package_import(module_name, display_name=None, required=True):
    """
    Test if a package can be imported
    
    Args:
        module_name: Name of the module to import
        display_name: Display name for the package (defaults to module_name)
        required: Whether this package is required (affects return value)
    
    Returns:
        tuple: (success: bool, version: str, error: str)
    """
    if display_name is None:
        display_name = module_name
    
    try:
        module = __import__(module_name)
        version = getattr(module, "__version__", "unknown")
        print(f"  ✓ {display_name:30} v{version}")
        return (True, version, None)
    except ImportError as e:
        status = "REQUIRED" if required else "OPTIONAL"
        print(f"  ✗ {display_name:30} FAILED ({status}): {e}")
        return (False, None, str(e))
    except Exception as e:
        print(f"  ✗ {display_name:30} ERROR: {e}")
        return (False, None, str(e))

def test_core_packages():
    """Test core scientific computing packages"""
    print_section("Core Scientific Packages")
    
    results = []
    packages = [
        ("numpy", "NumPy", True),
        ("pandas", "Pandas", True),
        ("scipy", "SciPy", True),
        ("sklearn", "scikit-learn", True),
    ]
    
    for pkg in packages:
        results.append((pkg[1], *test_package_import(*pkg)))
    
    return results

def test_machine_learning_packages():
    """Test machine learning packages"""
    print_section("Machine Learning Packages")
    
    results = []
    packages = [
        ("torch", "PyTorch", True),
        ("tensorflow", "TensorFlow", False),  # Optional
        ("xgboost", "XGBoost", True),
        ("lightgbm", "LightGBM", True),
    ]
    
    for pkg in packages:
        results.append((pkg[1], *test_package_import(*pkg)))
    
    return results

def test_fairness_packages():
    """Test fairness and bias mitigation packages"""
    print_section("Fairness & Bias Mitigation")
    
    results = []
    packages = [
        ("aif360", "AIF360", True),
    ]
    
    for pkg in packages:
        results.append((pkg[1], *test_package_import(*pkg)))
    
    return results

def test_dimensionality_reduction():
    """Test dimensionality reduction packages"""
    print_section("Dimensionality Reduction")
    
    results = []
    packages = [
        ("umap", "UMAP", True),
        ("opentsne", "openTSNE", True),
    ]
    
    for pkg in packages:
        results.append((pkg[1], *test_package_import(*pkg)))
    
    return results

def test_visualization_packages():
    """Test visualization packages"""
    print_section("Visualization Packages")
    
    results = []
    packages = [
        ("matplotlib", "Matplotlib", True),
        ("seaborn", "Seaborn", True),
        ("plotly", "Plotly", False),  # Optional
    ]
    
    for pkg in packages:
        results.append((pkg[1], *test_package_import(*pkg)))
    
    return results

def test_data_io_packages():
    """Test data I/O packages"""
    print_section("Data I/O Packages")
    
    results = []
    packages = [
        ("tables", "PyTables (HDF5)", True),
        ("h5py", "h5py", True),
        ("gdown", "gdown", True),
        ("osfclient", "OSF Client", True),
        ("pybiomart", "pyBioMart", True),
    ]
    
    for pkg in packages:
        results.append((pkg[1], *test_package_import(*pkg)))
    
    return results

def test_workflow_packages():
    """Test workflow management packages"""
    print_section("Workflow Management")
    
    results = []
    packages = [
        ("snakemake", "Snakemake", True),
    ]
    
    for pkg in packages:
        results.append((pkg[1], *test_package_import(*pkg)))
    
    return results

def test_utility_packages():
    """Test utility packages"""
    print_section("Utility Packages")
    
    results = []
    packages = [
        ("tqdm", "tqdm", True),
        ("psutil", "psutil", True),
        ("rich", "Rich", True),
        ("memory_profiler", "Memory Profiler", True),
        ("cpuinfo", "py-cpuinfo", True),
        ("jaxtyping", "jaxtyping", True),
        ("beartype", "beartype", True),
        ("tabulate", "tabulate", True),
    ]
    
    for pkg in packages:
        results.append((pkg[1], *test_package_import(*pkg)))
    
    return results

def test_numpy_functionality():
    """Test basic NumPy operations"""
    print_section("NumPy Functionality Test")
    
    try:
        import numpy as np
        
        # Create array
        arr = np.random.randn(1000, 100)
        print(f"  ✓ Created random array: {arr.shape}")
        
        # Basic operations
        mean = np.mean(arr)
        std = np.std(arr)
        print(f"  ✓ Statistics: mean={mean:.4f}, std={std:.4f}")
        
        # Linear algebra
        cov = np.cov(arr.T)
        print(f"  ✓ Covariance matrix: {cov.shape}")
        
        # SVD
        u, s, vt = np.linalg.svd(arr, full_matrices=False)
        print(f"  ✓ SVD decomposition: U{u.shape}, S{s.shape}, Vt{vt.shape}")
        
        return True
    except Exception as e:
        print(f"  ✗ NumPy functionality test failed: {e}")
        return False

def test_sklearn_functionality():
    """Test scikit-learn functionality"""
    print_section("scikit-learn Functionality Test")
    
    try:
        from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        from sklearn.metrics import accuracy_score, roc_auc_score
        import numpy as np
        
        # Create synthetic data
        np.random.seed(42)
        X = np.random.randn(500, 50)
        y = (X[:, 0] + X[:, 1] > 0).astype(int)
        
        # Test preprocessing
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        print(f"  ✓ StandardScaler: scaled data shape {X_scaled.shape}")
        
        # Test Random Forest
        rf = RandomForestClassifier(n_estimators=50, random_state=42, n_jobs=2)
        rf.fit(X_scaled, y)
        pred = rf.predict(X_scaled)
        acc = accuracy_score(y, pred)
        print(f"  ✓ Random Forest: accuracy={acc:.4f}")
        
        # Test Histogram Gradient Boosting
        hgb = HistGradientBoostingClassifier(max_iter=50, random_state=42)
        hgb.fit(X_scaled, y)
        pred_proba = hgb.predict_proba(X_scaled)[:, 1]
        auc = roc_auc_score(y, pred_proba)
        print(f"  ✓ HistGradientBoosting: AUC={auc:.4f}")
        
        # Test PCA
        pca = PCA(n_components=10)
        X_reduced = pca.fit_transform(X_scaled)
        var_explained = pca.explained_variance_ratio_.sum()
        print(f"  ✓ PCA: reduced to {X_reduced.shape[1]} dims, variance explained={var_explained:.4f}")
        
        return True
    except Exception as e:
        print(f"  ✗ scikit-learn functionality test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_torch_functionality():
    """Test PyTorch functionality"""
    print_section("PyTorch Functionality Test")
    
    try:
        import torch
        import torch.nn as nn
        
        # Create tensor
        x = torch.randn(100, 20)
        print(f"  ✓ Created tensor: {x.shape}")
        
        # Basic operations
        y = torch.matmul(x, x.T)
        print(f"  ✓ Matrix multiplication: {y.shape}")
        
        # Check CUDA availability
        cuda_available = torch.cuda.is_available()
        device_name = torch.cuda.get_device_name(0) if cuda_available else "CPU only"
        print(f"  ✓ Device: {device_name}")
        
        # Simple neural network
        model = nn.Sequential(
            nn.Linear(20, 10),
            nn.ReLU(),
            nn.Linear(10, 2)
        )
        output = model(x)
        print(f"  ✓ Neural network output: {output.shape}")
        
        # Test gradient computation
        loss = output.sum()
        loss.backward()
        print(f"  ✓ Gradient computation successful")
        
        return True
    except Exception as e:
        print(f"  ✗ PyTorch functionality test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_pandas_functionality():
    """Test Pandas functionality"""
    print_section("Pandas Functionality Test")
    
    try:
        import pandas as pd
        import numpy as np
        import tempfile
        
        # Create DataFrame
        df = pd.DataFrame({
            'A': np.random.randn(1000),
            'B': np.random.randn(1000),
            'C': np.random.choice(['X', 'Y', 'Z'], 1000),
            'D': np.random.randint(0, 100, 1000)
        })
        print(f"  ✓ Created DataFrame: {df.shape}")
        
        # Basic operations
        grouped = df.groupby('C').agg({'A': 'mean', 'B': 'std', 'D': 'sum'})
        print(f"  ✓ Grouped aggregation: {grouped.shape}")
        
        # Test CSV I/O
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            temp_csv = f.name
        
        try:
            df.to_csv(temp_csv, index=False)
            df_loaded = pd.read_csv(temp_csv)
            assert df_loaded.shape == df.shape
            print(f"  ✓ CSV I/O: wrote and read {df_loaded.shape}")
        finally:
            os.unlink(temp_csv)
        
        # Test HDF5 I/O
        with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as f:
            temp_h5 = f.name
        
        try:
            df.to_hdf(temp_h5, key='data', mode='w')
            df_loaded = pd.read_hdf(temp_h5, key='data')
            assert df_loaded.shape == df.shape
            print(f"  ✓ HDF5 I/O: wrote and read {df_loaded.shape}")
        finally:
            os.unlink(temp_h5)
        
        return True
    except Exception as e:
        print(f"  ✗ Pandas functionality test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_local_imports():
    """Test that local pipeline scripts can be imported"""
    print_section("Local Pipeline Imports")
    
    # Add scripts directory to path
    scripts_dir = Path(__file__).parent.parent / "scripts"
    if not scripts_dir.exists():
        print(f"  ⚠ Scripts directory not found at {scripts_dir}")
        return False
    
    sys.path.insert(0, str(scripts_dir))
    
    try:
        # Test utils module
        from utils import DataFrameCache, HashCache
        print(f"  ✓ utils.DataFrameCache")
        print(f"  ✓ utils.HashCache")
        
        # Test evaluations module
        from evaluations.util import repeated_cross_val
        print(f"  ✓ evaluations.util.repeated_cross_val")
        
        # Test prepdata module
        from prepdata.config import DATASETS_DIR
        print(f"  ✓ prepdata.config.DATASETS_DIR")
        
        return True
    except ImportError as e:
        print(f"  ✗ Failed to import local modules: {e}")
        import traceback
        traceback.print_exc()
        return False

def print_system_info():
    """Print system information"""
    print_header("System Information")
    
    print(f"  Hostname:          {platform.node()}")
    print(f"  Platform:          {platform.platform()}")
    print(f"  Python version:    {sys.version.split()[0]}")
    print(f"  Python executable: {sys.executable}")
    print(f"  Working directory: {os.getcwd()}")
    
    # Check if running on compute node
    hostname = platform.node()
    is_compute_node = any(x in hostname.lower() for x in ['compute', 'node', 'cn'])
    node_type = "Compute Node" if is_compute_node else "Login Node"
    print(f"  Node type:         {node_type}")
    
    # Environment variables
    print(f"\n  Environment Variables:")
    for var in ['PYTHON_ENV', 'UV_CACHE_DIR', 'SLURM_JOB_ID', 'SLURM_JOB_NAME']:
        value = os.environ.get(var, 'not set')
        print(f"    {var:20} = {value}")

def generate_summary(all_results, functionality_results):
    """Generate test summary"""
    print_header("Test Summary")
    
    # Count results
    total_packages = len(all_results)
    passed_packages = sum(1 for _, success, _, _ in all_results if success)
    failed_packages = total_packages - passed_packages
    
    total_functionality = len(functionality_results)
    passed_functionality = sum(1 for success in functionality_results if success)
    failed_functionality = total_functionality - passed_functionality
    
    print(f"\n  Package Imports:")
    print(f"    Total:   {total_packages}")
    print(f"    Passed:  {passed_packages} ({100*passed_packages/total_packages:.1f}%)")
    print(f"    Failed:  {failed_packages}")
    
    print(f"\n  Functionality Tests:")
    print(f"    Total:   {total_functionality}")
    print(f"    Passed:  {passed_functionality} ({100*passed_functionality/total_functionality:.1f}%)")
    print(f"    Failed:  {failed_functionality}")
    
    # List failures
    if failed_packages > 0:
        print(f"\n  Failed Package Imports:")
        for name, success, version, error in all_results:
            if not success:
                print(f"    ✗ {name}: {error}")
    
    # Overall status
    print()
    if failed_packages == 0 and failed_functionality == 0:
        print("  ✓ ALL TESTS PASSED")
        return 0
    else:
        print("  ✗ SOME TESTS FAILED")
        return 1

def main():
    """Run all validation tests"""
    start_time = datetime.now()
    
    print_header("Python Environment Validation")
    print(f"  Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Print system info
    print_system_info()
    
    # Test package imports
    print_header("Package Import Tests")
    
    all_results = []
    all_results.extend(test_core_packages())
    all_results.extend(test_machine_learning_packages())
    all_results.extend(test_fairness_packages())
    all_results.extend(test_dimensionality_reduction())
    all_results.extend(test_visualization_packages())
    all_results.extend(test_data_io_packages())
    all_results.extend(test_workflow_packages())
    all_results.extend(test_utility_packages())
    
    # Test functionality
    print_header("Functionality Tests")
    
    functionality_results = [
        test_numpy_functionality(),
        test_sklearn_functionality(),
        test_torch_functionality(),
        test_pandas_functionality(),
        test_local_imports(),
    ]
    
    # Generate summary
    exit_code = generate_summary(all_results, functionality_results)
    
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    
    print(f"\n  Completed: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Duration:  {duration:.2f} seconds")
    print("=" * 70)
    
    return exit_code

if __name__ == "__main__":
    sys.exit(main())
