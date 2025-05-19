import argparse
import numpy as np
import os
from util import *
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input-dir", type=Path, help="Input directory", required=True)
parser.add_argument("-o", "--output-path", type=Path, help="Path to output file", required=True)
args = parser.parse_args()

def calculate_mse(df1, df2):
    _, genes1 = split_discrete_continuous(df1)
    _, genes2 = split_discrete_continuous(df2)

    squared_error = (np.array(genes1) - np.array(genes2))**2

    return squared_error.mean()

# Open files
cache = DataFrameCache()

if not os.path.exists(args.output_path):
    with open(args.output_path, "w") as output_file:
        output_file.write("metric,adjuster,dataset,value\n")

# Calculate metrics
unadjusted_path = args.input_dir / "unadjusted.csv"
unadj = cache.get_dataframe(unadjusted_path)
dataset = os.path.basename(args.input_dir)

for method in ["scaled", "combat", "confounded"]:
    adjusted_path = args.input_dir / f"{method}.csv"
    df = cache.get_dataframe(adjusted_path)

    print(f"Calculating MSE for the {method} method on the {dataset} dataset", flush=True)
    value = calculate_mse(unadj, df)

    print(f"Saving output to {args.output_path}.")
    with open(args.output_path, "a") as output_file:
        output_file.write("{},{},{},{}\n".format("MSE", method, dataset, value))
