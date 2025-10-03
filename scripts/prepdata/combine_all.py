import pandas as pd
from pathlib import Path
import argparse

def print_now(*args, **kwargs):
    """Prints immediately to stdout (e.g., for real-time CLI feedback)."""
    print(*args, flush=True, **kwargs)

def combine_gold_unadjusted_files(gold_dir: Path, output_file: Path):
    """
    Combines all 'unadjusted.csv' files from GSE subfolders in 'gold_dir'.
    Each file gets a 'meta_source' column from its parent folder (e.g., GSE1234).
    Saves combined file to 'output_file'.
    """
    unadjusted_files = list(gold_dir.glob("gse*/unadjusted.csv"))
    if not unadjusted_files:
        raise FileNotFoundError(f"No unadjusted.csv files found in {gold_dir}/GSE*/")

    dfs = []
    for f in unadjusted_files:
        df = pd.read_csv(f, low_memory=False)
        gse_id = f.parent.name
        df['meta_source'] = gse_id
        print_now(f"Loaded {f} with shape {df.shape}")
        dfs.append(df)

    # Compute common columns across all files
    common_cols = set(dfs[0].columns)
    for df in dfs[1:]:
        common_cols.intersection_update(df.columns)
    common_cols = list(common_cols)

    print_now(f"Using {len(common_cols)} common columns across datasets.")

    combined_df = pd.concat([df[common_cols] for df in dfs], ignore_index=True)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    combined_df.to_csv(output_file, index=False)
    print_now(f"Saved combined data to {output_file} with shape {combined_df.shape}")

def main():
    parser = argparse.ArgumentParser(description="Combine all unadjusted.csv files from gold directory.")
    parser.add_argument('--input-dir', type=Path, required=True, help='Path to gold directory containing GSE subfolders')
    parser.add_argument('--output-file', type=Path, required=True, help='Path to save the combined CSV output')
    args = parser.parse_args()

    combine_gold_unadjusted_files(args.input_dir, args.output_file)

if __name__ == "__main__":
    main()