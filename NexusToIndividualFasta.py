#!/usr/bin/env python

from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import IUPAC,Gapped
import sys

#This script takes a nexus alignment and converts is to individual fasta files
#for each genome in the alignment

# check for correct arguments
if len(sys.argv) != 2:
    print("Usage: NexusToIndividualFasta.py <inputfile>")
    sys.exit(0)

input_name = sys.argv[1]

alignment = AlignIO.read(input_name, 'nexus')

for seq in alignment:
    out = open(seq.id + '.fasta', 'w')
    SeqIO.write(seq, out, 'fasta')
    out.close()

