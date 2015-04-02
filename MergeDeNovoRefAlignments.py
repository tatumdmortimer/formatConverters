#!/usr/bin/env python

from Bio import AlignIO
from Bio.Alphabet import IUPAC,Gapped
import sys, os, argparse

# This script takes an alignment of de novo assemblies to a reference and merges
# it with reference guided assemblies to the same reference sequence

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Combine whole genome \
with reference guided assembly alignment')
    parser.add_argument("wga", help="Fasta whole genome alignment", 
        action=FullPaths, type=is_file)
    parser.add_argument("rga", help="Nexus reference guided assembly alignment",
        action=FullPaths, type=is_file)
    parser.add_argument("reference", help="Name of reference in wga")
    parser.add_argument("outfile", help="Prefix for output alignment")
    return parser.parse_args()

def removeGaps(align, seqID):
    """ Removes positions in an alignment that are gaps in the provided \
    sequence"""
    refIndex = 0
    for seqRecord in align:
        if seqRecord.id == seqID:
            break
        refIndex += 1
    noGapList = []
    gapRegion = False
    noGapStart = 0
    noGapStop = 0
    for i in range(1, align.get_alignment_length()):
        if align[refIndex].seq[i] == '-':
            if not gapRegion:
                noGapStop = i
                noGapList.append((noGapStart, noGapStop))
                gapRegion = True
            else:
                continue
        elif i == align.get_alignment_length() - 1:
            noGapStop = i + 1
            noGapList.append((noGapStart, noGapStop))
        else:
            if gapRegion:
                noGapStart = i 
                gapRegion = False
            else:
                continue
    gapFreeAlign = align[:, noGapList[0][0]:noGapList[0][1]]
    for i in range(1,len(noGapList)):
        gapFreeAlign = gapFreeAlign + align[:, noGapList[i][0]:noGapList[i][1]]
    return gapFreeAlign 

args = get_args()

deNovoAlign = AlignIO.read(args.wga, 'fasta', alphabet=Gapped(IUPAC.ambiguous_dna, '-'))
refAlign = AlignIO.read(args.rga, 'nexus')

gapFreeAlign = removeGaps(deNovoAlign, args.reference)
print gapFreeAlign.get_alignment_length()

# add sequences from gap free de novo assembly alignment to reference guided
# assembly alignment (except for reference sequence)
for seq in gapFreeAlign:
    if seq.id != args.reference:
        refAlign.append(seq)

AlignIO.write(refAlign, args.outfile + ".nexus", 'nexus')
AlignIO.write(refAlign, args.outfile + ".fasta", 'fasta')

