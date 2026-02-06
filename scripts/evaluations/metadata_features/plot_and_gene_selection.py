#!/usr/bin/env python3
import os
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description="Select top genes per target from permutation importance across adjusters"
    )
    parser.add_argument("--importance_dir", required=True,
                        help="Directory with permutation importance CSVs per adjuster")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--top_k", type=int, default=20, help="Top k genes per target")
    parser.add_argument("--threshold", type=float, default=None,
                        help="Optional minimum importance (applied after top-k)")
    parser.add_argument("--agg", default="mean", choices=["mean", "max", "median"],
                        help="How to aggregate importance across adjusters")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # -------------------------
    # Load all CSVs
    # -------------------------
    all_files = [os.path.join(args.importance_dir, f)
                 for f in os.listdir(args.importance_dir) if f.endswith(".csv")]
    if not all_files:
        raise ValueError(f"No CSV files found in {args.importance_dir}")

    dfs = [pd.read_csv(f, index_col=0) for f in all_files]
    combined_df = pd.concat(dfs, axis=0)

    # Aggregate by gene
    if args.agg == "mean":
        agg_df = combined_df.groupby(combined_df.index).mean()
    elif args.agg == "max":
        agg_df = combined_df.groupby(combined_df.index).max()
    elif args.agg == "median":
        agg_df = combined_df.groupby(combined_df.index).median()

    selected_features = set()

    for target in agg_df.columns:
        ranked = agg_df[target].sort_values(ascending=False).head(args.top_k)
        if args.threshold is not None:
            ranked = ranked[ranked > args.threshold]

        selected_features.update(ranked.index)

        # Save per-target CSV
        per_target_csv = os.path.join(args.outdir, f"top_{args.top_k}_{target}_genes.csv")
        ranked.reset_index().rename(columns={"index": "gene", target: "importance"}).to_csv(per_target_csv, index=False)
        print(f"Saved top genes for {target}: {per_target_csv}")

    # Save union of selected genes
    filtered_df = agg_df.loc[sorted(selected_features)]
    filtered_csv = os.path.join(args.outdir, "filtered_permutation_importance.csv")
    filtered_df.to_csv(filtered_csv)
    print(f"Saved filtered permutation importance table: {filtered_csv}")

if __name__ == "__main__":
    main()
