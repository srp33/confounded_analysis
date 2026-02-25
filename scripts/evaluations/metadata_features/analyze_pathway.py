#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import gseapy as gp

def get_permutation_importance():

    pass

def run_preranked_gsea(ranked_file, output_dir, gmt_files):
    """
    Run preranked GSEA using a full .rnk file (tab-separated, no header).
    
    ranked_file : path to .rnk file (gene<TAB>score)
    output_dir  : folder to save GSEA results
    gmt_files   : list of local GMT file paths
    """
    # Load ranked file
    df = pd.read_csv(ranked_file, sep="\t", header=None, names=["gene", "score"])
    df["gene"] = df["gene"].str.upper()
    df = df.set_index("gene")
    ranking_metric = df["score"]

    for gmt_path in gmt_files:
        print(f"Running GSEA on {gmt_path}")
        gp.prerank(
            rnk=ranking_metric,
            gene_sets=gmt_path,
            outdir=os.path.join(output_dir, os.path.basename(gmt_path).replace(".gmt","")),
            permutation_num=100,  # adjust for speed/accuracy
            seed=42,
            verbose=True,
            min_size=5,
            max_size=1000
        )
        print(f"Completed {gmt_path}")

    print("All GSEA runs completed.")


def run_ora(selected_genes, gmt_files, outdir):
    """
    Run Over-Representation Analysis (ORA) using selected genes
    """
    res = gp.enrichr(
        gene_list=selected_genes,
        gene_sets=gmt_files,
        organism="Human",
        outdir=os.path.join(outdir, "ORA"),
        cutoff=0.05
    )

    if res.results is None or res.results.empty:
        print("No significant ORA terms found at cutoff=0.05")
        return pd.DataFrame()
    
    return res.results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full_ranked", required=True, help="Full ranked gene file for preranked GSEA (all genes)")
    parser.add_argument("--selected_genes_csv", required=True, help="Thresholded gene list for ORA (top genes)")
    parser.add_argument("--gmt_files", nargs="+", required=True, help="Paths to local GMT files")
    parser.add_argument("--outdir", required=True, help="Directory to save results")
    parser.add_argument("--top_n", type=int, default=None, help="Top N genes to use (optional)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load selected genes
    df = pd.read_csv(args.selected_genes_csv)

    df["gene"] = (
        df["gene"]
        .str.split(".").str[0]
        .str.upper()
    )

    if args.top_n is not None:
        df = df.nlargest(args.top_n, "max_importance")

    top_genes = df["gene"].tolist()

    print(f"Using {len(top_genes)} genes for pathway analysis")

    # Run GSEA
    run_preranked_gsea(
        ranked_file=args.full_ranked,
        output_dir=args.outdir,
        gmt_files=args.gmt_files
    )

    # Run ORA
    ora_results = run_ora(
        selected_genes=top_genes,
        gmt_files=args.gmt_files,
        outdir=args.outdir
    )

    # Save ORA results

    if not ora_results.empty:
        ora_results.to_csv(os.path.join(args.outdir, "ORA_results.csv"), index=False)
        print("ORA results saved.")
    else:
        print("No ORA results to save.")


if __name__ == "__main__":
    main()
