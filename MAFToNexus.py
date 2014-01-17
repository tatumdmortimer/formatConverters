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
    print("Usage: MAFToNexus.py <reference> <inputfile> <outputfile>")
    sys.exit(0)

ref = sys.argv[1]
input_name = sys.argv[2]
output_name = sys.argv[3]

input_file = open(input_name, 'r')
output_file = open(output_name, 'w')

# read in MAF file
multiple_alignment = AlignIO.parse(input_file, 'maf', alphabet=Gapped(IUPAC.ambiguous_dna, '-'))

# take blocks that contain reference sequence and add them to a list
alignment = []
for a in multiple_alignment:
    if a[0].id.startswith(ref):
        print a[0].annotations["start"]
        print int(a[0].annotations["start"]) + int(a[0].annotations["size"])
        alignment.append((int(a[0].annotations["start"]),a))

# sort the block list by start of ref sequence
sorted_alignment = sorted(alignment, key=itemgetter(0))

AlignIO.write(sorted_alignment[1][1], output_file, 'nexus')
input_file.close()
output_file.close()
