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
    parser.add_argument("--csvs", nargs="+", required=True,
                        help="List of permutation importance CSVs to aggregate")
    parser.add_argument("--outdir", required=True, help="Output directory for aggregated csv and heatmap")
    parser.add_argument("--agg", default="mean", choices=["mean", "max", "median"],
                        help="How to aggregate importance across adjusters")
    parser.add_argument("--target", required=True, help="Target name (used for filename)")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

        # -------------------------
    # Load all CSVs
    # -------------------------
    dfs = []
    for f in args.csvs:
        if not os.path.exists(f):
            raise FileNotFoundError(f"CSV not found: {f}")
        df = pd.read_csv(f, index_col=0)
        dfs.append(df)

    # Concatenate and aggregate
    combined_df = pd.concat(dfs, axis=0)
    if args.agg == "mean":
        agg_df = combined_df.groupby(combined_df.index).mean()
    elif args.agg == "median":
        agg_df = combined_df.groupby(combined_df.index).median()
    elif args.agg == "max":
        agg_df = combined_df.groupby(combined_df.index).max()

    # -------------------------
    # Save aggregated CSV
    # -------------------------
    agg_csv_path = os.path.join(args.outdir, f"{args.target}_aggregated_permutation_importance.csv")
    agg_df.to_csv(agg_csv_path)
    print(f"Saved aggregated CSV: {agg_csv_path}")

    # -------------------------
    # Plot heatmap
    # -------------------------
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        agg_df,
        cmap="viridis",
        annot=True,
        fmt=".3f",
        cbar_kws={"label": "Δ ROC AUC"}
    )
    plt.title(f"Permutation Importance Heatmap: {args.target}")
    plt.tight_layout()

    heatmap_path = os.path.join(args.outdir, f"{args.target}_permutation_importance_heatmap.png")
    plt.savefig(heatmap_path)
    plt.close()
    print(f"Saved heatmap: {heatmap_path}")

if __name__ == "__main__":
    main()