#!/usr/bin/env python

#import required libraries
import sys
sys.path.insert(1, "/home/peplab/src/alignio-maf")
try:
    from Bio import AlignIO
except ImportError:
    print "oops, the import didn't work"
from Bio.Alphabet import IUPAC,Gapped
from operator import itemgetter


#This script takes a MAF alignment and converts is to a
#nexus alignment

# check for correct arguments
if len(sys.argv) != 3:
    print("Usage: MAFToNexus.py <reference> <reference genome size> <inputfile> <outputfile>")
    sys.exit(0)

ref = sys.argv[1]
genomeSize = int(sys.argv[2])
input_name = sys.argv[3]
output_name = sys.argv[4]

infile = open(input_name, 'r')
outfile = open(output_name, 'w')

# read in MAF file
multiple_alignment = AlignIO.parse(infile, 'maf', alphabet=Gapped(IUPAC.ambiguous_dna, '-'))

# take blocks that contain reference sequence and add them to a list
alignment = []
for a in multiple_alignment:
    if a[0].id.startswith(ref):
        alignment.append((int(a[0].annotations["start"]),a))

# sort the block list by start of ref sequence
sorted_alignment = sorted(alignment, key=itemgetter(0))

# make a dictionary for genome sequences
genomeDict = {}

# read blocks and edit genome dictionary
for block in sorted_alignment:
    block = block[1] # get rid of tuple format    
    # get start and stop sites for reference
    blockStart = int(a[0].annotations["start"])
    blockStop = blockStart + int(a[0].annotations["size"])
    for seqRecord in block:
        strainName = seqRecord.id.split('.')[0]
        if strainName not in genomeDict:
            # add sequence to genome dictionary
            # using '-' as placeholder until actual sequence is added
            genomeDict[strainName] = ['-']*genomeSize
        genome = genomeDict[strainName]
        genome[blockStart:blockStop] = list(seqRecord.seq)
        genomeDict[strainName] = genome
        
outfile.close()
