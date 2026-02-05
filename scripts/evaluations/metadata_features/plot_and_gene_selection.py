# Analysis and visualization

# Load permutation importance CSV
# Filtering logic

import os 
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -------------------------
# Main
# -------------------------
def main():
    parser = argparse.ArgumentParser(description="Plot permutation importance and select top genes")
    parser.add_argument("--importance_csv", required=True, help="Permutation importance CSV from File 2")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--top_k", type=int, default=20, help="Select top k features per target")
    parser.add_argument("--threshold", type=float, default=0.005, help="Minimum importance to keep feature")
    
    args = parser.parse_args()

    # -------------------------
    # Load permutation importance
    # -------------------------
    imp_df = pd.read_csv(args.importance_csv, index_col=0)
    
    # -------------------------
    # Filter features (keep if above threshold in any target)
    # -------------------------
    filtered_df = imp_df[(imp_df > args.threshold).any(axis=1)]
    
    # -------------------------
    # Save filtered table
    # -------------------------
    os.makedirs(args.outdir, exist_ok=True)
    filtered_csv = os.path.join(args.outdir, "filtered_permutation_importance.csv")
    filtered_df.to_csv(filtered_csv)
    print(f"Saved filtered permutation importance: {filtered_csv}")

    # -------------------------
    # Generate heatmap
    # -------------------------
    plt.figure(figsize=(10, 8))
    sns.heatmap(filtered_df, annot=True, cmap="viridis", fmt=".3f")
    plt.title("Filtered Permutation Importance per Target")
    
    heatmap_path = os.path.join(args.outdir, "permutation_importance_heatmap.png")
    plt.tight_layout()
    plt.savefig(heatmap_path)
    plt.close()
    print(f"Saved heatmap: {heatmap_path}")

    # -------------------------
    # Select top K features per target
    # -------------------------
    top_features_dict = {}
    for target in filtered_df.columns:
        top_genes = filtered_df[target].sort_values(ascending=False).head(args.top_k)
        top_features_dict[target] = top_genes.index.tolist()

        # Save per-target CSV
        per_target_csv = os.path.join(args.outdir, f"top_{args.top_k}_{target}_genes.csv")
        pd.DataFrame({
            "gene": top_genes.index,
            "importance": top_genes.values
        }).to_csv(per_target_csv, index=False)
        print(f"Saved top {args.top_k} genes for {target}: {per_target_csv}")

    print("Gene selection complete.")

if __name__ == "__main__":
    main()