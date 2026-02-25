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
    parser.add_argument("--csvs", nargs="+", required=True,
        help="List of permutation importance CSVs (one per adjuster)"
    )
    parser.add_argument("--outdir",required=True,
        help="Output directory for heatmaps"
    )
    parser.add_argument( "--threshold",type=float,default=0.005,
        help="Filter genes with importance > threshold in at least one adjuster"
    )
    parser.add_argument("--test",required=True,
        help="Name of the test set"
    )
    parser.add_argument("--train",required=True,
        help="Name of the train set"
    )
    parser.add_argument("--aligned_csv",required=True,
        help="CSV file of already aligned metadata for the training and testing sets"
    )
    parser.add_argument("--adjusted_csvs",nargs="+",required=True,
        help="List of CSVs of the adjusted data"
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
    print(f"Targets: {targets}")
    
    # -------------------------
    # Sanity check: compare target columns
    # -------------------------
    print("\n=== Column diagnostics ===")

    # Collect column sets
    column_sets = {adj: set(df.columns) for adj, df in dfs.items()}

    # Print columns per adjuster
    for adj, cols in column_sets.items():
        print(f"\nAdjuster: {adj}")
        print(f"  Number of targets: {len(cols)}")
        print(f"  Targets: {sorted(cols)}")

    # Use first adjuster as reference
    reference_adj = next(iter(column_sets))
    reference_cols = column_sets[reference_adj]

    print(f"\nReference adjuster: {reference_adj}")

    # Compare each adjuster to reference
    for adj, cols in column_sets.items():
        if adj == reference_adj:
            continue

        missing_in_adj = reference_cols - cols
        extra_in_adj = cols - reference_cols

        if missing_in_adj or extra_in_adj:
            print(f"\nColumn mismatch in adjuster: {adj}")

            if missing_in_adj:
                print(f"  Missing columns ({len(missing_in_adj)}):")
                print(f"    {sorted(missing_in_adj)}")

            if extra_in_adj:
                print(f"  Extra columns ({len(extra_in_adj)}):")
                print(f"    {sorted(extra_in_adj)}")

    # If you still want strict enforcement:
    all_columns_equal = all(cols == reference_cols for cols in column_sets.values())
    if not all_columns_equal:
        raise ValueError("Column mismatch detected across adjusters.")

    # -------------------------
    # Plot one heatmap per target
    # -------------------------

    # Collect genes that pass threshold with their maximum importance
    gene_max_importance = {}
        
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
        
        # Update max importance for each gene
        for gene, row in target_df.iterrows():
            max_val = row.max()
            if gene in gene_max_importance:
                gene_max_importance[gene] = max(gene_max_importance[gene], max_val)
            else:
                gene_max_importance[gene] = max_val

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

    # -------------------------
    # Save selected genes (thresholded) and full ranked file (all genes)
    # -------------------------

    # Thresholded genes for ORA
    selected_genes = list(gene_max_importance.keys())

    if selected_genes:
        print(f"\nSubsetting Datasets to {len(selected_genes)} genes...")

        # Save thresholded dataset subsets
        subset_dataset(
            input_csv=args.aligned_csv,
            selected_genes=selected_genes,
            outdir=args.outdir,
            train_name=args.train,
            test_name=args.test,
            prefix="unadjusted"
        )

        for adjusted_csv in args.adjusted_csvs:
            subset_dataset(
                input_csv=adjusted_csv,
                selected_genes=selected_genes,
                outdir=args.outdir,
                train_name=args.train,
                test_name=args.test,
                prefix="adjusted"
            )

        # Save thresholded gene list
        gene_list_df = pd.DataFrame({
            "gene": selected_genes,
            "max_importance": [gene_max_importance[g] for g in selected_genes]
        }).sort_values("max_importance", ascending=False)

        gene_list_path = os.path.join(args.outdir, "selected_genes.csv")
        gene_list_df.to_csv(gene_list_path, index=False)
        print(f"Saved thresholded gene list for ORA: {gene_list_path}")

    # -------------------------
    # Save full ranked gene list for preranked GSEA
    # -------------------------

    def compute_max_gene_importance(dfs):
        """
        Compute maximum gene importance across all adjusters and all targets.
        Returns a pd.Series sorted descending.
        """
        combined = pd.concat(dfs.values(), axis=1)  # genes × all targets from all adjusters
        gene_max = combined.max(axis=1).sort_values(ascending=False)
        return gene_max

    gene_rank_series = compute_max_gene_importance(dfs)

    gene_rank_df = pd.DataFrame({
        "gene": gene_rank_series.index,
        "score": gene_rank_series.values
    })

    rank_path = os.path.join(args.outdir, "gsea_prerank_max_importance.rnk")

    # Tab-separated, no header, suitable for gseapy prerank
    gene_rank_df.to_csv(rank_path, sep="\t", index=False, header=False)
    print(f"Saved full ranked GSEA file (all genes): {rank_path}")
        
def subset_dataset(input_csv, selected_genes, outdir, train_name, test_name, prefix):
    print(f"\nProcessing dataset: {input_csv}")

    df = pd.read_csv(input_csv, low_memory=False)

    if "meta_source" not in df.columns:
        raise ValueError("meta_source column not found in dataset.")

    # Identify metadata columns
    meta_cols = [col for col in df.columns if col.startswith("meta_")]

    print(f"Keeping {len(meta_cols)} metadata columns:")
    print(meta_cols)

    # Keep only genes that exist in dataset
    available_genes = [g for g in selected_genes if g in df.columns]

    missing_genes = set(selected_genes) - set(available_genes)
    if missing_genes:
        print(f"Warning: {len(missing_genes)} selected genes not found in dataset.")

    cols_to_keep = meta_cols + available_genes
    df_subset = df[cols_to_keep]

    train_df = df_subset[df_subset["meta_source"] == train_name]
    test_df = df_subset[df_subset["meta_source"] == test_name]
    
    print(f"Train samples: {train_df.shape[0]}")
    print(f"Test samples: {test_df.shape[0]}")

    base = os.path.splitext(os.path.basename(input_csv))[0]

    train_out = os.path.join(outdir, f"{prefix}_{base}_train_selected_genes.csv")
    test_out = os.path.join(outdir, f"{prefix}_{base}_test_selected_genes.csv")

    train_df.to_csv(train_out, index=False)
    test_df.to_csv(test_out, index=False)

    print(f"Saved train subset: {train_out}")
    print(f"Saved test subset: {test_out}")

if __name__ == "__main__":
    main()