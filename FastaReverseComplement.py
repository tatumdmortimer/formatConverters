#!/usr/bin/env python

from Bio import AlignIO
from Bio.Alphabet import IUPAC,Gapped
import sys

#This script takes a FASTA alignment and converts is to a
#nexus alignment

# check for correct arguments
if len(sys.argv) != 3:
    print("Usage: FastaReverseComplement.py <inputfile> <outputfile>")
    sys.exit(0)

input_name = sys.argv[1]
output_name = sys.argv[2]


records = [rec.reverse_complement(id="rc_"+rec.id, description = "reverse
complement") for rec in SeqIO.parse(input_name, "fasta")]
SeqIO.write(records, output_name, "fasta")
AlignIO.write(alignment, output_file, 'nexus')

