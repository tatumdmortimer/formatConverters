#!/usr/bin/env python

from Bio import AlignIO
import sys

#This script takes a sequential phylip alignment and converts is to a
#FASTA alignment


input_name = sys.argv[1]
output_name = sys.argv[2]

input_file = open(input_name, 'r')
output_file = open(output_name, 'w')

alignment = AlignIO.read(input_file, 'phylip-sequential')
AlignIO.write(alignment, output_file, 'phylip')

input_file.close()
output_file.close()
