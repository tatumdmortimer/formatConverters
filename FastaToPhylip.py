#!/usr/bin/env python

from Bio import AlignIO
from Bio.Alphabet import IUPAC,Gapped
import sys

#This script takes a FASTA alignment and converts is to a
#phylip sequential alignment

# check for correct arguments
if len(sys.argv) != 3:
    print("Usage: FastaToPhylip.py <inputfile> <outputfile>")
    sys.exit(0)

input_name = sys.argv[1]
output_name = sys.argv[2]

input_file = open(input_name, 'r')
output_file = open(output_name, 'w')

alignment = AlignIO.read(input_file, 'fasta', alphabet=Gapped(IUPAC.ambiguous_dna, '-'))
AlignIO.write(alignment, output_file, 'phylip-sequential')

input_file.close()
output_file.close()
