import pandas as pd
import os
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input-path", type=Path, help="Path to original log file (csv)", required=True)
parser.add_argument("-o", "--output-path", type=Path, help="Path to output file", required=True)
args = parser.parse_args()

# Input: csv, where each row is: metric,adjuster,dataset,value\n
# Output: Columns are adjuster names: scale, combat, etc. Rows are dataset names. Cells are values.

pivoted = pd.read_csv(args.input_path)\
            .drop("metric", axis=1)\
            .pivot(index='dataset', columns='adjuster', values='value')\
            .round(3)
   
pivoted.to_csv(args.output_path)