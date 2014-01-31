#!/usr/bin/env python

from Bio import AlignIO
from Bio.Alphabet import IUPAC,Gapped
import sys

# This script takes an alignment of de novo assemblies to a reference and merges
# it with reference guided assemblies to the same reference sequence

def removeGaps(align, seqID):
    """ Removes positions in an alignment that are gaps in the provided \
    sequence"""
    refIndex = 0
    for seqRecord in align:
        if seqRecord.id == seqID:
            break
        refIndex += 1
    gapFreeAlign = align[:, 0:1]
    for i in range(1, align.get_alignment_length()):
        if align[refIndex].Seq[0] != '-':
            gapFreeAlign = gapFreeAlign + align[:, i:i+1]
    return gapFreeAlign 

# check for correct arguments
if len(sys.argv) != 5:
    print("Usage: FastaToNexus.py <denovo file> <refguided file> <reference \
    name> <outputfile>")
    sys.exit(0)

deNovoName = sys.argv[1]
refFileName = sys.argv[2]
reference = sys.argv[3]
outFileName = sys.argv[4]

deNovoFile = open(deNovoName, 'r')
refFile = open(refFileName, 'r')
outFile = open(outFileName. 'w')

deNovoAlign = AlignIO.read(deNovoFile, 'fasta', 
    alphabet=Gapped(IUPAC.ambiguous_dna, '-'))
refAlign = AlignIO.read(refFile, 'nexus')

gapFreeAlign = removeGaps(deNovoAlign, reference)

# add sequences from gap free de novo assembly alignment to reference guided
# assembly alignment (except for reference sequence)
for seq in gapFreeAlign:
    if seq.id != reference:
        refAlign.append(seq)

AlignIO.write(refAlign, outFile, 'nexus')

deNovoFile.close()
refFile.close()
outFile.close()
