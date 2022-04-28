import argparse
import os
import numpy as np
import sys
import time
output = "learningRateTestResults.csv"
if not os.path.exists(output):
    with open(output, "w") as output_file:
        output_file.write("sampleSize, numRandomFeatures, numRedundantFeatures, numInformativeFeatures, adjuster, codeSize, batchResult, trueResult, sclaer, weightLoss, learningRate\n")
        #output_file.write("Result\n")

#file = open(sys.argv[1], "r")
#values = file.split()
#print(file.readline()) #results
results = []
fileNum = 1
learningRate = ['0.1', '1.0', '5.0', '10']
weightLoss = ['1.0']
scaler = ['linear']
codesamples = ['100']
samples = ['100', '200', '300']
randomFeatures = ['800', '950', '990']
redundantFeatures = ['175', '40', '8']
informativeFeatures = ['25', '10', '2']
adjusters = ['unadjusted','scaled','combat','confounded']

for sampleSize in samples:
    for codeSize in codesamples:
        for scaler in scalers:
            for weight in weightLoss:
                for place in range(3):

                    batch = open(sys.argv[1] + "s" + sampleSize + "+f" + informativeFeatures[place] + "+c" + codeSize + "+" + scaler + "+w" + weight + ".csv", "r")
                    true = open(sys.argv[2] + "s" + sampleSize + "+f" + informativeFeatures[place] + "+c" + codeSize + "+" + scaler + "+w" + weight +  ".csv", "r")
            batchvalues = []
            truevalues = []

            for line in batch:
                cells = line.split(",")
                batchvalues.append(cells[3].strip())
            for line in true:
                cells = line.split(",")
                truevalues.append(cells[3].strip())

            for adjust in range(4):
                batchResult = batchvalues[adjust + 2]
                trueResult = truevalues[adjust + 2]
                codeNum = codeSize
                results.append([sampleSize, randomFeatures[place], redundantFeatures[place], informativeFeatures[place], adjusters[adjust], codeSize, batchResult, trueResult, scaler, weight])
with open(output, 'a') as output_file:
    for line in results:
        output_file.write(",".join(line) + "\n")

#file.close()
