import argparse
import os
import pandas as pd
from sklearn.svm import SVC
import os.path
from os import path
import numpy as np

import sys
import time
import random

INPUT_DIR = "../data/"
output = "singlemetriccomparison_minus.csv"

def random_accuracy(data, label) :
  samples  = open(INPUT_DIR + data + "/unadjusted.csv")
  out = []
  for line in samples :
      cells = line.split(",")
      out.append(cells)
  index = [i for i,x in enumerate(out[0]) if x == label]
  #print(index)
  out = np.array(out)
  #out = out[1:,]
  out  = out[:,int(index[0])]
  #print(out)
  a, cnts = np.unique(out, return_counts=True)
  high_freq, high_freq_element = cnts.max(), a[cnts.argmax()]
  max = high_freq
  total = len(out) - 1
  #print(samples)
  return 1.0 * max / total


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input-dir", help="Input directory", required=True)
parser.add_argument("-o", "--output-dir", help="Path to output directory", required=True)
args = parser.parse_args()


if not os.path.exists(args.output_dir + output):
    with open(args.output_dir + output, "w") as output_file:
        output_file.write("metric, adjuster, dataset, score\n")
        #output_file.write("Result\n")

datasets = ["simulated_expression", "bladderbatch", "gse37199", 
                "tcga_small", "tcga_medium", "tcga"
                ]
batch_labels = ["Batch", "batch", "plate", 
              "CancerType", "CancerType", "CancerType" 
              ]
true_labels = ["Class", "batch", "Stage", 
                "TP53_Mutated", "TP53_Mutated", "TP53_Mutated"
                ]
results = []

if path.exists(args.input_dir + "batch_classification.csv"):

    batch_raw = open(args.input_dir + "batch_classification.csv")
    true_raw = open(args.input_dir + "true_classification.csv")
    batch = []
    true = []
    for line in batch_raw:
        line = line.strip()
        elements = line.split(',')
        batch.append(elements)
    for line in true_raw :
        line = line.strip()
        elements = line.split(',')
        true.append(elements)
    batch = np.array(batch)
    true = np.array(true)
    for i in range(len(batch)):
        if batch[i,2] == "bladderbatch" :
            true = np.insert(true, i, np.array([batch[i,0], batch[i,1], batch[i,2], '0']), axis = 0)
    selections = np.array([False, False, False, True])
    batch_values = batch[:, selections]
    all_data = np.hstack((true, batch_values))

    for ime in range(len(datasets)) :
        
        batch_random = random_accuracy(datasets[ime], batch_labels[ime])
        true_random = random_accuracy(datasets[ime], true_labels[ime])
        for line in all_data:
            if line[2] == datasets[ime] :
                #score = (float(line[3]) - true_random) + ( batch_random - float(line[4]))
                #score = float(line[3]) - abs(float(line[4]) - batch_random)
                score = (float(line[3]) - true_random) - abs(batch_random - float(line[4]))
                results.append([line[0], line[1], line[2], str(score)])


with open(args.output_dir + output, 'a') as output_file:
    for line in results:
        output_file.write(",".join(line) + "\n")


        
