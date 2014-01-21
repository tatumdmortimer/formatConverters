#!/usr/bin/env python

from Bio import AlignIO
from StringIO import StringIO
from Bio.Nexus import Nexus
from Bio.Alphabet import IUPAC,Gapped
import sys

#This script takes a FASTA alignment and converts is to a
#nexus alignment

# check for correct arguments
if len(sys.argv) != 3:
    print("Usage: FastaToNexusInterleaved.py <inputfile> <outputfile>")
    sys.exit(0)

input_name = sys.argv[1]
output_name = sys.argv[2]

alignment = AlignIO.read(input_name, 'fasta', alphabet=Gapped(IUPAC.ambiguous_dna, '-'))

output = StringIO()
AlignIO.write(alignment, output, 'nexus')
p = Nexus.Nexus()
p.read(output.getvalue())
p.write_nexus_data(output_name, interleave=True)

