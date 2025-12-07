# This script:
# 1. Reads the combined_csv
# 2. Extracts each unique GSE id
# 3. Produces one order file for the specified test source

import os
import pandas as pd
import random
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Generate order file for a specific test source from combined data."
    )
    parser.add_argument("--input", required=True, help="Path to combined CSV file")
    parser.add_argument("--output", required=True, help="Path to output order file")
    parser.add_argument("--test-source", required=False, help="Test source ID (optional, will be inferred from output filename)")
    
    args = parser.parse_args()
    
    all_combined_csv = args.input
    output_file = args.output
    output_dir = os.path.dirname(output_file)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")
        
    if not os.path.exists(all_combined_csv):
        print(f"File not found: {all_combined_csv}")
        raise FileNotFoundError(f"File not found: {all_combined_csv}")

    df = pd.read_csv(all_combined_csv)
ue()

    # Extract test source from argument or filename
    :
        test_source t_source
    
        test_source = os.path.basenam
    
    random.seed(234)
 
    # Split into training and testing
    train_source = [x for xe_ids if  gs in
        "train_source": train_source
    })

    # Write CSV
    out_df.to_csv(output_file, index=False)
    print(f"Wrote: {output_file}")

if __name__ == "__main__":
    main()
        

        