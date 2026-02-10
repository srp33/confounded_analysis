#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description="Plot per-target permutation importance heatmaps across batch adjusters"
    )
    parser.add_argument(
        "--csvs",
        nargs="+",
        required=True,
        help="List of permutation importance CSVs (one per adjuster)"
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for heatmaps"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.005,
        help="Filter genes with importance > threshold in at least one adjuster"
    )

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # -------------------------
    # Load CSVs (one per adjuster)
    # -------------------------
    dfs = {}

    for f in args.csvs:
        if not os.path.exists(f):
            raise FileNotFoundError(f"CSV not found: {f}")

        df = pd.read_csv(f, index_col=0)

        # infer adjuster name from directory
        adjuster = os.path.basename(os.path.dirname(f))
        dfs[adjuster] = df

    # -------------------------
    # Sanity check: same targets everywhere
    # -------------------------
    targets = dfs[next(iter(dfs))].columns.tolist()

    for adj, df in dfs.items():
        if list(df.columns) != targets:
            raise ValueError(f"Target mismatch in adjuster: {adj}")

    # -------------------------
    # Plot one heatmap per target
    # -------------------------
    for target in targets:

        # genes × adjusters
        target_df = pd.DataFrame({
            adjuster: dfs[adjuster][target]
            for adjuster in dfs
        })

        # Filter genes
        target_df = target_df[(target_df > args.threshold).any(axis=1)]

        if target_df.empty:
            print(f"Skipping {target}: no genes pass threshold")
            continue

        plt.figure(figsize=(10, max(6, 0.25 * target_df.shape[0])))

        sns.heatmap(
            target_df,
            cmap="viridis",
            annot=False,
            cbar_kws={"label": "Δ ROC AUC"},
            linewidths=0.2
        )

        plt.title(f"Permutation Importance – {target}")
        plt.xlabel("Batch adjustment method")
        plt.ylabel("Gene")
        plt.tight_layout()

        out_path = os.path.join(
            args.outdir,
            f"permutation_importance_{target}.png"
        )
        plt.savefig(out_path, dpi=300)
        plt.close()

        print(f"Saved heatmap: {out_path}")


if __name__ == "__main__":
    main()
