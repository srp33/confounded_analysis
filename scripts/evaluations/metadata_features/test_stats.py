import argparse
import os
from pathlib import Path

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, linregress


def run_tests(file_path, meta_cols, outdir):
    df = pd.read_csv(file_path)

    # Validate metadata columns
    for col in meta_cols:
        if col not in df.columns:
            raise ValueError(f"Metadata column '{col}' not found in {file_path}")

    metadata = df[meta_cols]

    # Identify gene columns
    gene_cols = [c for c in df.columns if c not in meta_cols and not c.startswith("meta")]
    genes = df[gene_cols].apply(pd.to_numeric, errors='coerce')

    # Output DataFrames
    tstat_df = pd.DataFrame(index=gene_cols)
    slope_df = pd.DataFrame(index=gene_cols)
    pval_df = pd.DataFrame(index=gene_cols)

    base_name = Path(file_path).stem
    os.makedirs(outdir, exist_ok=True)

    for meta_col in meta_cols:
        labels = metadata[meta_col].dropna()

        # Align genes with non-null metadata
        valid_idx = labels.index
        gene_data = genes.loc[valid_idx]
        labels = labels.loc[valid_idx]

        unique_vals = labels.unique()

        # -------------------------
        # Binary → Welch t-test
        # -------------------------
        if len(unique_vals) == 2:
            group1 = gene_data[labels == unique_vals[0]]
            group2 = gene_data[labels == unique_vals[1]]

            if group1.shape[0] < 2 or group2.shape[0] < 2:
                print(f"[Skipping t-test] {meta_col} (not enough samples)")
                continue

            t_stats, p_vals = ttest_ind(
                group1,
                group2,
                axis=0,
                equal_var=False,
                nan_policy="omit"
            )

            tstat_df[meta_col] = t_stats

        # -------------------------
        # Continuous → Linear regression
        # -------------------------
        else:
            slopes = []
            pvals = []

            x = labels.values.astype(float)

            for gene in gene_cols:
                y = gene_data[gene].values

                # Remove NaNs pairwise
                mask = ~np.isnan(x) & ~np.isnan(y)
                if mask.sum() < 3:
                    slopes.append(np.nan)
                    pvals.append(np.nan)
                    continue

                try:
                    res = linregress(x[mask], y[mask])
                    slopes.append(res.slope)
                    pvals.append(res.pvalue)
                except Exception:
                    slopes.append(np.nan)
                    pvals.append(np.nan)

            slope_df[meta_col] = slopes
            pval_df[meta_col] = pvals

    # -------------------------
    # Save outputs
    # -------------------------
    tstat_path = os.path.join(outdir, f"{base_name}-tstats.csv")
    slope_path = os.path.join(outdir, f"{base_name}-slopes.csv")
    pval_path = os.path.join(outdir, f"{base_name}-pvalues.csv")

    if not tstat_df.empty:
        tstat_df.index.name = "Gene"
        tstat_df.to_csv(tstat_path)
        print(f"[Saved] {tstat_path}")

    if not slope_df.empty:
        slope_df.index.name = "Gene"
        slope_df.to_csv(slope_path)
        print(f"[Saved] {slope_path}")

    if not pval_df.empty:
        pval_df.index.name = "Gene"
        pval_df.to_csv(pval_path)
        print(f"[Saved] {pval_path}")


def main():
    parser = argparse.ArgumentParser(description="Gene-wise tests for metadata columns.")
    parser.add_argument("--csv", required=True, help="Input CSV file")
    parser.add_argument("--meta-cols", nargs="+", required=True, help="Metadata columns to test")
    parser.add_argument("--outdir", required=True, help="Output directory")

    args = parser.parse_args()

    print(f"\nProcessing {args.csv}...")
    run_tests(args.csv, args.meta_cols, args.outdir)


if __name__ == "__main__":
    main()