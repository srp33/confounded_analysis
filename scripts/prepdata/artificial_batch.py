import argparse
import pandas as pd
import numpy as np
from random import seed
from utils.adjustments import Noise, split_discrete_continuous

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Path to input file")
parser.add_argument("outfile", help="Path to output file")
args = parser.parse_args()

seed(0)
df = pd.read_csv(args.infile)
discrete, continuous = split_discrete_continuous(df)
n_pixels = len(continuous.columns)
n1 = Noise((n_pixels,), order=10, discount_factor=0.02)
n2 = Noise((n_pixels,), order=10, discount_factor=0.02)
batch = []
images = []
for i, row in continuous.iterrows():
    add_noise = np.random.uniform(-0.1, 0.1, size=(n_pixels,))
    mult_noise = np.random.uniform(0.85, 0.95, size=(n_pixels,))
    image = row.tolist()
    image = (image + add_noise) * mult_noise
    if i % 2 == 0:
        image = n1.adjust(image)
        batch.append("A")
    else:
        image = n2.adjust(image)
        batch.append("B")
    images.append(image)
noisy = pd.DataFrame(images)
stuff = pd.concat([discrete, noisy], axis="columns")
cols = stuff.columns.tolist()
stuff["Batch"] = batch
stuff = stuff[["Batch"] + cols]
stuff.to_csv(args.outfile, index=False)
