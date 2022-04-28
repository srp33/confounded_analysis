import argparse
import os
import numpy as np
import sys
import time
output = "confoundedTestResults.csv"
if not os.path.exists(output):
    with open(output, "w") as output_file:
        output_file.write("numSamples, numRandomFeatures, numRedundantFeatures, numInformativeFeatures, adjuster, batchResult, trueResult, TestNumber\n")
        #output_file.write("Result\n")

#file = open(sys.argv[1], "r")
#values = file.split()
#print(file.readline()) #results
results = []
fileNum = 1
samples = ['100','200','300']
randomFeatures = ['800','950','990']
redundantFeatures = ['175','40','8']
informativeFeatures = ['25','10','2']
adjusters = ['unadjusted','scaled','combat','confounded']

for sample in samples:
    for feature in range(3):
        batch = open(sys.argv[1] + str(fileNum) + ".csv", "r")
        true = open(sys.argv[2] + str(fileNum) + ".csv", "r")
        batchvalues = []
        truevalues = []

        for line in batch:
            cells = line.split(",")
            batchvalues.append(cells[3].strip())
        for line in true:
            cells = line.split(",")
            truevalues.append(cells[3].strip())

        for adjust in range(4):
            numSamples = sample
            numRandomFeatures = randomFeatures[feature]
            numRedundantFeatures = redundantFeatures[feature]
            numInformativeFeatures = informativeFeatures[feature]
            batchResult = batchvalues[adjust + 2]
            trueResult = truevalues[adjust + 2]
            results.append([numSamples, numRandomFeatures, numRedundantFeatures, numInformativeFeatures, adjusters[adjust], batchResult, trueResult, str(fileNum)])
        fileNum += 1
with open(output, 'a') as output_file:
    for line in results:
        output_file.write(",".join(line) + "\n")

#file.close()