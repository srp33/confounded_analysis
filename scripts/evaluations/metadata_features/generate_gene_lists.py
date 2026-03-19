import pandas as pd
import numpy as np
import argparse
import os

# -------------------------
# Load CSVs
# -------------------------
def load_csvs(paths):
    dfs = []
    for p in paths:
        df = pd.read_csv(p, index_col=0)
        df = df.apply(pd.to_numeric, errors="coerce")
        dfs.append(df)
    return dfs

# -------------------------
# Align genes across adjusters
# -------------------------
def intersect_genes(dfs):
    common = set(dfs[0].index)
    for df in dfs[1:]:
        common &= set(df.index)
    common = sorted(common)

    if len(common) == 0:
        raise ValueError("No common genes across adjusters.")

    print(f"Genes after intersection: {len(common)}")
    return [df.loc[common] for df in dfs], common

# -------------------------
# Aggregate across adjusters
# -------------------------
def aggregate(dfs, method="max_abs"):
    stacked = np.stack([df.values for df in dfs])

    if method == "max_abs":
        agg_vals = np.nanmax(np.abs(stacked), axis=0)

    elif method == "mean":
        agg_vals = np.nanmean(stacked, axis=0)

    elif method == "median":
        agg_vals = np.nanmedian(stacked, axis=0)

    else:
        raise ValueError(f"Unknown aggregation method: {method}")

    return pd.DataFrame(
        agg_vals,
        index=dfs[0].index,
        columns=dfs[0].columns
    )

# -------------------------
# Transform scores (for regression / p-values)
# -------------------------
def transform_scores(df, mode="tstat"):
    if mode == "tstat":
        return df  # already appropriate

    elif mode == "slope":
        return df  # keep sign

    elif mode == "abs_slope":
        return df.abs()

    elif mode == "pval":
        # Convert to -log10(p)
        return -np.log10(df)

    else:
        raise ValueError(f"Unknown mode: {mode}")

# -------------------------
# Save top genes
# -------------------------
def save_top_genes(df, outdir, threshold=None):
    os.makedirs(outdir, exist_ok=True)

    for target in df.columns:
        ranked = df[target].sort_values(ascending=False)

        if threshold is not None:
            ranked = ranked[ranked > threshold]

        out_path = os.path.join(outdir, f"{target}_top_genes.csv")
        ranked.to_csv(out_path, header=[target])
        print(f"Saved top genes: {out_path}")

# -------------------------
# Save GSEA files
# -------------------------
def save_rnk_files(df, outdir):
    os.makedirs(outdir, exist_ok=True)

    for target in df.columns:
        ranked = df[target].sort_values(ascending=False)

        out_path = os.path.join(outdir, f"{target}.rnk")
        ranked.to_csv(out_path, sep="\t", header=False)
        print(f"Saved GSEA file: {out_path}")

# -------------------------
# Main
# -------------------------
def main(args):
    dfs = load_csvs(args.csvs)
    dfs, genes = intersect_genes(dfs)

    agg_df = aggregate(dfs, method=args.agg)

    score_df = transform_scores(agg_df, mode=args.mode)

    save_top_genes(score_df, args.outdir, threshold=args.threshold)
    save_rnk_files(score_df, args.outdir)

# -------------------------
# CLI
# -------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--csvs", nargs="+", required=True)
    parser.add_argument("--outdir", required=True)

    parser.add_argument(
        "--mode",
        choices=["tstat", "slope", "abs_slope", "pval"],
        default="tstat",
        help="Type of input data"
    )

    parser.add_argument(
        "--agg",
        choices=["max_abs", "mean", "median"],
        default="max_abs",
        help="Aggregation across adjusters"
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=None
    )

    args = parser.parse_args()
    main(args)