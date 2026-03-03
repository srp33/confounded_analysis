import argparse
import os
from pathlib import Path

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

def run_ttests(file_path, meta_cols, outdir):
    """
    Compute gene-wise Welch t-tests for multiple metadata columns.
    Save a single CSV per adjuster: genes x metadata columns.
    """
    df = pd.read_csv(file_path)

    # Check metadata columns exist
    for col in meta_cols:
        if col not in df.columns:
            raise ValueError(f"Metadata column '{col}' not found in {file_path}")

    # Separate metadata and genes
    metadata = df[meta_cols]

    # Exclude columns starting with 'meta' that are NOT in meta_cols
    gene_cols = [c for c in df.columns if c not in meta_cols and not c.startswith("meta")]
    genes = df[gene_cols].apply(pd.to_numeric, errors='coerce')  # ensure numeric

    # Prepare DataFrame to hold all t-statistics
    tstat_df = pd.DataFrame(index=gene_cols)

    base_name = Path(file_path).stem
    os.makedirs(outdir, exist_ok=True)

    for meta_col in meta_cols:
        labels = metadata[meta_col]
        unique_vals = labels.dropna().unique()

        if len(unique_vals) != 2:
            print(f"[Skipping] {meta_col} (not binary)")
            continue

        group1 = genes[labels == unique_vals[0]]
        group2 = genes[labels == unique_vals[1]]

        if group1.shape[0] < 2 or group2.shape[0] < 2:
            print(f"[Skipping] {meta_col} (not enough samples)")
            continue

        # Compute t-statistics
        t_stats, _ = ttest_ind(
            group1,
            group2,
            axis=0,
            equal_var=False,
            nan_policy="omit"
        )

        tstat_df[meta_col] = t_stats

    # Save one CSV per adjuster
    out_file = os.path.join(outdir, f"{base_name}_tstats.csv")
    tstat_df.index.name = "Gene"
    tstat_df.to_csv(out_file)
    print(f"[Saved] {out_file}")


def main():
    parser = argparse.ArgumentParser(description="Compute gene-wise t-tests for multiple metadata columns.")
    parser.add_argument("files", nargs="+", help="Input CSV files (samples as rows)")
    parser.add_argument("--meta-cols", nargs="+", required=True, help="Metadata columns to test")
    parser.add_argument("--outdir", required=True, help="Directory to store t_test results")
    args = parser.parse_args()

    for f in args.files:
        print(f"\nCalculating t-statistics for {f}...")
        run_ttests(f, args.meta_cols, args.outdir)

if __name__ == "__main__":
    main()