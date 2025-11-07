import pandas as pd
import glob
import os

results_root = "results"
for adjuster_dir in os.listdir(results_root):
    adjuster_path = os.path.join(results_root, adjuster_dir)
    if os.path.isdir(adjuster_path):
        all_files = glob.glob(os.path.join(adjuster_path, "*_metrics.csv"))
        df = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
        df.to_csv(os.path.join(results_root, f"{adjuster_dir}_metrics_aggregated.csv"), index=False)
        print(f"Aggregated results for {adjuster_dir} -> {adjuster_dir}_metrics_aggregated.csv")
