#!/usr/bin/env python3
"""
Test script for Python-only execution via run_with_env.sh
Tests that Python environment is activated correctly
"""

import sys
import os

def main():
    print("=== Python Environment Test ===")
    print(f"Python version: {sys.version}")
    print(f"Python executable: {sys.executable}")
    print(f"Virtual env: {os.environ.get('VIRTUAL_ENV', 'Not set')}")
    
    # Test that we can import key packages
    try:
        import numpy as np
        print(f"✓ NumPy {np.__version__} imported successfully")
    except ImportError as e:
        print(f"✗ NumPy import failed: {e}")
        return 1
    
    try:
        import pandas as pd
        print(f"✓ Pandas {pd.__version__} imported successfully")
    except ImportError as e:
        print(f"✗ Pandas import failed: {e}")
        return 1
    
    try:
        import sklearn
        print(f"✓ scikit-learn {sklearn.__version__} imported successfully")
    except ImportError as e:
        print(f"✗ scikit-learn import failed: {e}")
        return 1
    
    print("\n✓ All Python tests passed!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
