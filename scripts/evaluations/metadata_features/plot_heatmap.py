import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

def load_adjuster_data(paths):
    dfs = {}
    for path in paths:
        adjuster = os.path.basename(os.path.dirname(path))
        df = pd.read_csv(path, index_col=0)
        dfs[adjuster] = df
    return dfs

def plot_heatmaps(dfs, outdir):
    os.makedirs(outdir, exist_ok=True)

    # assume all adjusters share same target columns
    targets = list(next(iter(dfs.values())).columns)

    for target in targets:
        # genes × adjusters matrix
        matrix = pd.DataFrame({
            adj: dfs[adj][target]
            for adj in dfs
        })

        matrix = matrix.dropna()

        plt.figure(figsize=(10, max(6, 0.25 * matrix.shape[0])))

        sns.heatmap(
            matrix,
            cmap="viridis",
            cbar_kws={"label": target}
        )

        plt.title(f"Heatmap for {target}")
        plt.xlabel("Adjuster")
        plt.ylabel("Genes")
        plt.tight_layout()

        out_path = os.path.join(outdir, f"{target}_heatmap.png")
        plt.savefig(out_path, dpi=300)
        plt.close()

        print(f"Saved: {out_path}")