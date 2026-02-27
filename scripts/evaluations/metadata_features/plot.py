import argparse
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -------------------------
# Utilities
# -------------------------
def clean_genes(genes: pd.Index | list) -> pd.Index:
    """Normalize gene names: strip whitespace, uppercase."""
    return pd.Index([str(g).strip().upper() for g in genes])

def load_expression_matrix(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = clean_genes(df.index)
    return df

def load_ranked_genes(path: str, top_n: int) -> pd.Index:
    ranked = pd.read_csv(path, index_col=0)
    genes = ranked.index[:top_n]
    return clean_genes(genes)

def subset_expression(X: pd.DataFrame, genes: pd.Index) -> pd.DataFrame:
    genes_present = [g for g in genes if g in X.index]
    if len(genes_present) < len(genes):
        missing = set(genes) - set(genes_present)
        print(f"Warning: {len(missing)} genes not found in expression matrix: {list(missing)[:5]} ...")
    return X.loc[genes_present]

def plot_heatmap(df: pd.DataFrame, title: str, out_path: str):
    plt.figure(figsize=(8, max(4, len(df)//5)))
    sns.heatmap(df, cmap="vlag", center=0)
    plt.title(title)
    plt.ylabel("Genes")
    plt.xlabel("Samples")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
    print(f"Saved heatmap: {out_path}")

# -------------------------
# Main
# -------------------------
def main(args):
    os.makedirs(args.outdir, exist_ok=True)

    expr_df = load_expression_matrix(args.expression_csv)

    for ranked_path in args.ranked_lists:
        ranked_name = os.path.basename(ranked_path).replace("_ranked.csv", "")
        genes = load_ranked_genes(ranked_path, top_n=args.top_n)

        df_subset = subset_expression(expr_df, genes)
        if df_subset.empty:
            print(f"Skipping {ranked_name}: no genes found in expression matrix.")
            continue

        out_path = os.path.join(args.outdir, f"{ranked_name}_heatmap.png")
        plot_heatmap(df_subset, title=f"Top {args.top_n} genes – {ranked_name}", out_path=out_path)

# -------------------------
# CLI
# -------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression_csv", required=True)
    parser.add_argument("--ranked_lists", nargs="+", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--top_n", type=int, default=50)
    args = parser.parse_args()

    main(args)