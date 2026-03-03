import pandas as pd
import numpy as np
import argparse
import os

def load_csvs(paths):
    return [pd.read_csv(p, index_col=0) for p in paths]

def compute_max_across_adjusters(dfs):
    """
    Restrict to genes present in ALL adjusters (intersection),
    then compute max absolute importance across adjusters.
    """
    common_genes = set(dfs[0].index)
    for df in dfs[1:]:
        common_genes &= set(df.index)
    common_genes = sorted(common_genes)

    print(f"Number of genes after intersection: {len(common_genes)}")

    if len(common_genes) == 0:
        raise ValueError("No common genes across adjusters.")
    
    aligned = [df.loc[common_genes] for df in dfs]
    stacked = np.stack([df.values for df in aligned])
    max_vals = np.max(np.abs(stacked), axis=0)

    return pd.DataFrame(
        max_vals,
        index=common_genes,
        columns=aligned[0].columns,
    )

def save_top_genes(df, outdir, threshold=None):
    """
    Save top genes above threshold for each target
    """
    os.makedirs(outdir, exist_ok=True)

    for target in df.columns:
        ranked = df[target].abs().sort_values(ascending=False)
        if threshold is not None:
            ranked = ranked[ranked > threshold]
        out_path = os.path.join(outdir, f"{target}_top_genes_ttest.csv")
        ranked.to_csv(out_path, header=[target])
        print(f"Top genes saved for {target}: {out_path}")

def save_rnk_files(df, outdir):
    """
    Save one .rnk file per target for GSEA (gene + score)
    """
    os.makedirs(outdir, exist_ok=True)

    for target in df.columns:
        ranked = df[target].sort_values(ascending=False)
        out_path = os.path.join(outdir, f"{target}_ttest.rnk")
        ranked.to_csv(out_path, sep="\t", header=False)
        print(f"Ranked file for GSEA saved: {out_path}")

def main(args):
    dfs = load_csvs(args.csvs)
    max_df = compute_max_across_adjusters(dfs)

    # Save top genes per target (above threshold)
    save_top_genes(max_df, args.outdir, threshold=args.threshold)

    # Save one .rnk file per target for GSEA
    save_rnk_files(max_df, args.outdir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--csvs", nargs="+", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--threshold", type=float, default=None)
    args = parser.parse_args()

    main(args)