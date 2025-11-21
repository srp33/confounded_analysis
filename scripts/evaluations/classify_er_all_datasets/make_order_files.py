# This R/Python script should:
# 1. Read the combined_csv
# 2. Extract each unique GSE id
# 3. Produce one file per GSE id in gse_order_files/

import os
import sys
import pandas as pd
import random

def main():
    if len(sys.argv) != 3: 
        print("Usage: python make_order_files.py <combined_csv> <output_dir>")
        sys.exit(1)

    all_combined_csv = sys.argv[1]
    output_dir = sys.argv[2]

    if not os.path.exists(all_combined_csv):
        print(f"File not found: {all_combined_csv}")
        sys.exit(1)

    df = pd.read_csv(all_combined_csv)

    gse_ids = df['meta_source'].unique()

    random.seed(234)

    for id in gse_ids: 
        # Split into training and testing
        test_source = id
        train_source = [x for x in gse_ids if x != test_source]
        random.shuffle(train_source)

        # Build a dataframe for output
        out_df = pd.DataFrame({
            "train_source": train_source
        })

        out_file = os.path.join(output_dir, f"{test_source}_order.csv")

        # Write CSV
        out_df.to_csv(out_file, index=False)

        print(f"Wrote: {out_file}")

if __name__ == "__main__":
    main()
        

        