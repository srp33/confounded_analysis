#!/usr/bin/env python3
import os
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description="Select top genes per target from permutation importance CSVs"
    )
    parser.add_argument(
        "--csvs",
        nargs="+",
        required=True,
        help="List of permutation importance CSVs to aggregate"
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--top_k", type=int, default=20, help="Top k genes per target")
    parser.add_argument("--threshold", type=float, default=None, help="Minimum importance")
    parser.add_argument(
        "--agg",
        default="mean",
        choices=["mean", "median", "max"],
        help="Aggregation across adjusters"
    )

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

    combined_df = pd.concat(dfs, axis=0)

    # Aggregate by gene
    if args.agg == "mean":
        agg_df = combined_df.groupby(combined_df.index).mean()
    elif args.agg == "median":
        agg_df = combined_df.groupby(combined_df.index).median()
    elif args.agg == "max":
        agg_df = combined_df.groupby(combined_df.index).max()

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
