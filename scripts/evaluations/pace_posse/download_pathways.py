#!/usr/bin/env python3
"""
Download pathways locally to cache them for POSSE
"""

import sys
import os
sys.path.insert(0, 'scripts')

import gseapy as gp
import pickle

def download_pathways():
    """Download and cache pathway databases"""
    
    print("Downloading MSigDB Hallmark pathways...")
    
    try:
        # Download the pathways
        pathways = gp.get_library(name='MSigDB_Hallmark_2020', organism='Human')
        
        print(f"Successfully downloaded {len(pathways)} pathways")
        
        # Save to cache file
        cache_file = 'scripts/pathways_cache.pkl'
        with open(cache_file, 'wb') as f:
            pickle.dump(pathways, f)
        
        print(f"Pathways cached to {cache_file}")
        
        # Print some pathway info
        print("\nSample pathways:")
        for i, (name, genes) in enumerate(list(pathways.items())[:5]):
            print(f"  {name}: {len(genes)} genes")
        
        return True
        
    except Exception as e:
        print(f"Error downloading pathways: {e}")
        return False

if __name__ == "__main__":
    success = download_pathways()
    if success:
        print("✅ Pathways downloaded successfully!")
    else:
        print("❌ Failed to download pathways")
        sys.exit(1)