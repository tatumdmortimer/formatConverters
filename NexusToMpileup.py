#!/usr/bin/env python

#This script will generate an mpileup formatted file from a "nexus" alignment file of sequences.

import sys
from collections import defaultdict

def usage():
    print "script.py <nexusAlignment> <outFileName>"

dataDict = defaultdict(list)
def make_dataDict(nexusAlignment):
    #make a dictionary with position as key and character as value
    with open(nexusAlignment, 'r') as inputDat:
        for line in inputDat:
            line = line.strip()
            if len(line) < 100:
               continue 
            else:
                seq = list(line)
                for i, c in enumerate(seq):
                    if c != '-':
                        position = i + 1
                        dataDict[position].append(c)
    return dataDict

def write_outfile(dataDict, outFileName):
    with open(outFileName, "w") as outfile:
        for key in range(4411533):
            if key in dataDict and len(dataDict[key]) > 1:
                ref = dataDict[key][0]
                alleles = dataDict[key][1:]
                numStrains = len(alleles)
                conf = [5] * numStrains
                outfile.write("%s\t%i\t%s\t%i\t%s\t%s\n" %
                ("gi|57116681|ref|NC_000962.2|",
                key,
                ref,
                numStrains,
                "".join(alleles),
                "".join(str(i) for i in conf))
                )         

if len(sys.argv) < 2:
    usage()
    sys.exit(0)

nexusAlignment, outFileName = sys.argv[1:]
dataDict = make_dataDict(nexusAlignment)
write_outfile(dataDict, outFileName)
