#!/usr/bin/env python3
import argparse
import pandas as pd
import os

def load_ranked(path):
    df = pd.read_csv(path, index_col=0)
    return df

def main(args):
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    all_genes = {}

    for path in args.ranked_lists:
        df = load_ranked(path)
        col = df.columns[0]

        df = df.sort_values(by=col, ascending=False)

        if args.top_n:
            df = df.head(args.top_n)

        if args.threshold is not None:
            df = df[df[col] > args.threshold]

        for gene, val in df[col].items():
            if gene in all_genes:
                all_genes[gene] = max(all_genes[gene], val)
            else:
                all_genes[gene] = val

    out_df = pd.DataFrame({
        "gene": list(all_genes.keys()),
        "score": list(all_genes.values())
    }).sort_values("score", ascending=False)

    out_df.to_csv(args.output, index=False)
    print(f"Saved gene list: {args.output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranked_lists", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--top_n", type=int, default=None)
    parser.add_argument("--threshold", type=float, default=None)
    args = parser.parse_args()
    main(args)