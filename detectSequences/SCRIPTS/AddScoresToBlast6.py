# Fix tabular format
import sys
import os
import time

stime = time.time()

score = None
scoreDict = {}
rows = []
for n, line in enumerate(open(sys.argv[1])):
    if n<=1: continue
    line = line.strip().split()
    if len(line) > 0:
        if line[0] =="a":
            score = float(line[1].replace("score=", ""))
            rows = []
            continue
        elif line[0] == "s":
            rows.append(line[0:6])
            continue
    elif len(line) == 0 and score is not None:
            try:
                evalue = (float(rows[0][5]) * float(rows[1][5])) / (2 ** score)
            except:
                evalue = 0.0
            if rows[1][4] == "-":
                fixBp = int(rows[1][5]) - 1 - int(rows[1][2])
                scoreDict["{}_{}_{}_{}".format(rows[0][1], int(rows[0][2]), rows[1][1], fixBp )] = (evalue, score)
            else:
                scoreDict["{}_{}_{}_{}".format(rows[0][1], int(rows[0][2]), rows[1][1], int(rows[1][2]))] = (evalue, score)
            score = ""
    else: 
        continue

for line in open(sys.argv[2]):
    line = line.strip().split()
    evalue, score = scoreDict["{}_{}_{}_{}".format(line[1], int(line[8])-1, line[0], int(line[6]) - 1)]
    if score > 0: print("{}\t{:.2e}\t{}".format('\t'.join(line), evalue, int(score)))