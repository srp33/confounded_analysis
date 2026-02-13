# Run pre-ranked GSEA
# gene -> max/mean permutation importance across adjusters

# Run ORA on special genes
# Genes stable across >= 3 adjusters

# Compare pathway enrichment consistency across adjusters

#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import numpy as np
import gseapy as gp


def compute_max_importance(df):
    return df.max(axis=1)


def run_preranked_gsea(ranking_series, adjuster_name, gene_sets, outdir):
    print(f"Running pre-ranked GSEA for {adjuster_name}")

    rnk_df = (
        ranking_series
        .sort_values(ascending=False)
        .reset_index()
    )
    rnk_df.columns = ["gene", "score"]

    rnk_path = os.path.join(outdir, f"{adjuster_name}_prerank.rnk")
    rnk_df.to_csv(rnk_path, sep="\t", index=False, header=False)

    res = gp.prerank(
        rnk=rnk_path,
        gene_sets=gene_sets,
        permutation_num=1000,
        outdir=os.path.join(outdir, f"{adjuster_name}_gsea"),
        seed=42,
        verbose=True
    )

    return res.res2d

def run_ora(selected_genes, background_genes, adjuster_name, gene_sets, outdir):
    """Run ORA (over-representation analysis)."""
    print(f"Running ORA for {adjuster_name}")

    res = gp.enrichr(
        gene_list=selected_genes,
        gene_sets=gene_sets,
        background=background_genes,
        organism="Human",
        outdir=os.path.join(outdir, f"{adjuster_name}_ora"),
        cutoff=0.05
    )

    return res.results


def compare_pathway_consistency(results_dict, fdr_threshold=0.05):
    """
    Compare pathway enrichment consistency across adjusters.
    Returns a pathway × adjuster matrix of significant enrichments.
    """
    print("Comparing pathway consistency...")

    all_pathways = set()
    for df in results_dict.values():
        all_pathways.update(df.index)

    consistency = pd.DataFrame(
        0,
        index=sorted(all_pathways),
        columns=results_dict.keys()
    )

    for adjuster, df in results_dict.items():
        sig = df[df["FDR q-val"] < fdr_threshold]
        consistency.loc[sig.index, adjuster] = 1

    return consistency


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--perm_csvs", nargs="+", required=True)
    parser.add_argument("--selected_genes_csv", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Choose pathway databases
    gene_sets = [
        "Hallmark_2020",
        "KEGG_2021_Human",
        "Reactome_2022"
    ]

    # Load selected genes (for ORA)
    selected_genes = pd.read_csv(args.selected_genes_csv)["gene"].tolist()

    prerank_results = {}
    ora_results = {}

    for csv in args.perm_csvs:

        adjuster = os.path.basename(os.path.dirname(csv))
        print(f"\nProcessing adjuster: {adjuster}")

        df = pd.read_csv(csv, index_col=0)

        # Full ranking for GSEA
        ranking = compute_max_importance(df)

        # Background genes
        background_genes = ranking.index.tolist()

        # ---- Run Pre-ranked GSEA ----
        gsea_df = run_preranked_gsea(
            ranking_series=ranking,
            adjuster_name=adjuster,
            gene_sets=gene_sets,
            outdir=args.outdir
        )

        prerank_results[adjuster] = gsea_df

        # ---- Run ORA ----
        ora_df = run_ora(
            selected_genes=selected_genes,
            background_genes=background_genes,
            adjuster_name=adjuster,
            gene_sets=gene_sets,
            outdir=args.outdir
        )

        ora_results[adjuster] = ora_df.set_index("Term")

    # -----------------------------------------------------
    # Compare Consistency (GSEA)
    # -----------------------------------------------------

    gsea_consistency = compare_pathway_consistency(prerank_results)
    gsea_consistency_path = os.path.join(args.outdir, "gsea_pathway_consistency.csv")
    gsea_consistency.to_csv(gsea_consistency_path)

    print(f"Saved GSEA consistency matrix: {gsea_consistency_path}")

    # -----------------------------------------------------
    # Compare Consistency (ORA)
    # -----------------------------------------------------

    ora_consistency = compare_pathway_consistency(ora_results)
    ora_consistency_path = os.path.join(args.outdir, "ora_pathway_consistency.csv")
    ora_consistency.to_csv(ora_consistency_path)

    print(f"Saved ORA consistency matrix: {ora_consistency_path}")


if __name__ == "__main__":
    main()