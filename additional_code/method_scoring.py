import argparse
import os
import numpy as np
import sys
import time

file = open(sys.argv[1], "r")

batchUnadjusted = 0
batchScaled = 0
batchCombat = 0
batchConfounded = 0

trueUnadjusted = 0
trueScaled = 0
trueCombat = 0
trueConfounded = 0

numBatch = 0
numTrue = 0
for line in file:
    cells = line.split(",")
    if (cells[4] == "unadjusted"):
        batchUnadjusted += float(cells[6])
        trueUnadjusted += float(cells[7])
        numBatch += 1
        numTrue += 1
    if (cells[4] == "scaled"):
        batchScaled += float(cells[6])
        trueScaled += float(cells[7])
    if (cells[4] == "combat"):
        batchCombat += float(cells[6])
        trueCombat += float(cells[7])
    if (cells[4] == "confounded"):
        batchConfounded += float(cells[6])
        trueConfounded += float(cells[7])
print("data from " + sys.argv[1] + "\n")
print("Method Scores for code Size as 100:\n")
print("Unadjusted\tbatch: " + str(batchUnadjusted/numBatch) + "\ttrue: " + str(trueUnadjusted/numTrue) + "\n")
print("Scaled\t\tbatch: " + str(batchScaled/numBatch) + "\ttrue: " + str(trueScaled/numTrue) + "\n")
print("Combat\t\tbatch: " + str(batchCombat/numBatch) + "\ttrue: " + str(trueCombat/numTrue) + "\n")
print("Confounded\tbatch: " + str(batchConfounded/numBatch) + "\ttrue: " + str(trueConfounded/numTrue) + "\n")
#with open(output, 'a') as output_file:
    #for line in results:
        #output_file.write(",".join(line) + "\n")

#file.close()