import os
import argparse
import pandas as pd
import gseapy as gp
from glob import glob

def run_preranked_gsea(ranked_file, output_dir, gmt_files):
    df = pd.read_csv(ranked_file, sep="\t", header=None, names=["gene", "score"])
    df["gene"] = df["gene"].str.upper()
    df = df.set_index("gene")
    ranking_metric = df["score"]

    for gmt_path in gmt_files:
        outdir_gmt = os.path.join(output_dir, os.path.basename(gmt_path).replace(".gmt",""))
        os.makedirs(outdir_gmt, exist_ok=True)
        print(f"Running GSEA on {gmt_path}")
        gp.prerank(
            rnk=ranking_metric,
            gene_sets=gmt_path,
            outdir=outdir_gmt,
            permutation_num=100,  # adjust for speed/accuracy
            seed=42,
            verbose=True,
            min_size=5,
            max_size=1000
        )
        print(f"Completed {gmt_path}")

def run_ora(selected_genes, gmt_files, outdir):
    """
    Run Over-Representation Analysis (ORA) using selected genes
    """
    if len(selected_genes) < 2:
        print(f"Skipping ORA: too few genes ({len(selected_genes)})")
        return pd.DataFrame()

    try:
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
    except ValueError as e:
        print(f"Skipping ORA due to error: {e}")
        return pd.DataFrame()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene_lists_dir", required=True, help="Directory containing per-target _top_genes.csv and .rnk files")
    parser.add_argument("--gmt_files", nargs="+", required=True, help="Paths to local GMT files")
    parser.add_argument("--outdir", required=True, help="Directory to save results")
    parser.add_argument("--top_n", type=int, default=None, help="Top N genes to use (optional)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Loop over all targets based on _top_genes.csv
    top_gene_files = glob(os.path.join(args.gene_lists_dir, "*_top_genes.csv"))

    for top_csv in sorted(top_gene_files):
        target_name = os.path.basename(top_csv).replace("_top_genes.csv", "")
        rnk_file = os.path.join(args.gene_lists_dir, f"{target_name}.rnk")

        if not os.path.exists(rnk_file):
            print(f"Skipping {target_name}: no .rnk file found")
            continue

        # Load top genes for ORA
        df_top = pd.read_csv(top_csv)
        df_top["gene"] = df_top.iloc[:, 0].str.upper()  # first column = gene names
        top_genes = df_top["gene"].tolist()
        if args.top_n is not None:
            top_genes = top_genes[:args.top_n]

        target_outdir = os.path.join(args.outdir, target_name)
        os.makedirs(target_outdir, exist_ok=True)

        print(f"\nRunning pathway analysis for target: {target_name}")
        print(f"Top genes: {len(top_genes)}, .rnk file: {rnk_file}")

        # Run GSEA
        run_preranked_gsea(rnk_file, target_outdir, args.gmt_files)

        # Run ORA
        ora_results = run_ora(top_genes, args.gmt_files, target_outdir)
        if not ora_results.empty:
            ora_results.to_csv(os.path.join(target_outdir, "ORA_results.csv"), index=False)
            print(f"ORA results saved for {target_name}")

    print("\nAll targets processed.")

if __name__ == "__main__":
    main()