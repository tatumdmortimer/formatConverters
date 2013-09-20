#!/usr/bin/env python

from Bio import AlignIO
import sys

#This script takes a sequential phylip alignment and converts is to a
#FASTA alignment

# check for correct arguments
if len(sys.argv) != 3:
    print("Usage: PhylipSequentialToFasta.py <inputfile> <outputfile>")
    sys.exit(0)


input_name = sys.argv[1]
output_name = sys.argv[2]

input_file = open(input_name, 'r')
output_file = open(output_name, 'w')

alignment = AlignIO.read(input_file, 'phylip-sequential')
AlignIO.write(alignment, output_file, 'fasta')

input_file.close()
output_file.close()
