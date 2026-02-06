#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(
        description="Plot heatmap of permutation importance from multiple adjuster CSVs"
    )
    parser.add_argument("--importance_dir", required=True,
                        help="Directory containing permutation importance CSVs for each adjuster")
    parser.add_argument("--outdir", required=True, help="Output directory for heatmap")
    parser.add_argument("--agg", default="mean", choices=["mean", "max", "median"],
                        help="How to aggregate importance across adjusters")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # -------------------------
    # Load all CSVs and concatenate
    # -------------------------
    all_files = [
        os.path.join(args.importance_dir, f) 
        for f in os.listdir(args.importance_dir) 
        if f.endswith(".csv")
    ]

    if not all_files:
        raise ValueError(f"No CSV files found in {args.importance_dir}")

    dfs = []
    for f in all_files:
        df = pd.read_csv(f, index_col=0)
        dfs.append(df)

    # Align on genes
    combined_df = pd.concat(dfs, axis=0)
    
    # Aggregate by gene (row-wise)
    if args.agg == "mean":
        agg_df = combined_df.groupby(combined_df.index).mean()
    elif args.agg == "max":
        agg_df = combined_df.groupby(combined_df.index).max()
    elif args.agg == "median":
        agg_df = combined_df.groupby(combined_df.index).median()

    # -------------------------
    # Heatmap
    # -------------------------
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        agg_df,
        cmap="viridis",
        annot=True,
        fmt=".3f",
        cbar_kws={"label": "Mean Δ ROC AUC"}
    )
    plt.title("Permutation Importance Across Adjusters")
    plt.tight_layout()

    heatmap_path = os.path.join(args.outdir, "permutation_importance_heatmap.png")
    plt.savefig(heatmap_path)
    plt.close()
    print(f"Saved heatmap: {heatmap_path}")

if __name__ == "__main__":
    main()
