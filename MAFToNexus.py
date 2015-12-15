#!/usr/bin/env python

#import required libraries
import sys
sys.path.insert(1, "/opt/PepPrograms/alignio-maf")
try:
    from Bio import AlignIO
    from Bio.AlignIO import MafIO
except ImportError:
    print "oops, the import didn't work"
from Bio.Alphabet import IUPAC,Gapped
from Bio.Align import MultipleSeqAlignment
from operator import itemgetter


#This script takes a MAF alignment and converts is to a
#nexus alignment

# check for correct arguments
if len(sys.argv) != 6:
    print("Usage: MAFToNexus.py <reference> <reference genome size> \
    <total number of sequences> <inputfile> <outputfile>")
    sys.exit(0)

reference = sys.argv[1]
genomeSize = int(sys.argv[2])
numSeq = int(sys.argv[3])
input_name = sys.argv[4]
output_prefix = sys.argv[5]

infile = open(input_name, 'r')
outfile = open(output_prefix, 'w')

# read in MAF file
multiple_alignment = list(AlignIO.parse(infile, 'maf', 
alphabet=Gapped(IUPAC.ambiguous_dna, '-')))
refList = []
refGenome = [0]*genomeSize
# get alignments with all sequences and reverse complement those that
# are on incorrect strand
for a in multiple_alignment:
    strainIDs = []
    for seqRecord in a:
        strainIDs.append(seqRecord.id)
    if len(strainIDs) == numSeq:
        refInd = strainIDs.index(reference)
        if a[refInd].seq.count('-') == len(a[refInd].seq):
            print "block"
            continue
        elif a[refInd].annotations["strand"] == "-1":
            start = int(a[refInd].annotations["srcSize"]) - int(a[refInd].annotations["start"]) - int(a[refInd].annotations["size"])
            stop = start + int(a[refInd].annotations["size"])
            for i in range(start,stop):
                refGenome[i] += 1
                if refGenome[i] > 1:
                    print i
            newAlignment = []
            for seqRecord in a:
                seq = seqRecord.reverse_complement(id=True, name=True,
                    features=True, annotations=True)
                if seq.annotations["strand"] == "-1":
                    seq.annotations["strand"] = "+1"
                    seq.annotations["start"] = (int(seq.annotations["srcSize"]) 
                    - int(seq.annotations["start"]) -
                    int(seq.annotations["size"]))
                elif seq.annotations["strand"] == "+1":
                    seq.annotations["strand"] = "-1"
                    seq.annotations["start"] = int(seq.annotations["srcSize"]) -int(seq.annotations["start"]) + int(seq.annotations["size"]) + 1
                newAlignment.append(seq)
            refList.append(MultipleSeqAlignment(newAlignment))
        else:
            start = int(a[refInd].annotations["start"])
            stop = start + int(a[refInd].annotations["size"])
            for i in range(start,stop):
                refGenome[i] += 1
                if refGenome[i] > 1:
                    print i
            refList.append(a)

print refGenome.count(0)
print refGenome.count(1)
print refGenome.count(2)
print refGenome.count(3)
print refGenome.count(4)

# create new temporary MAF file
AlignIO.write(refList, input_name + ".tmp", "maf")

# make MAF index file
idx = MafIO.MafIndex(reference + ".mafindex", input_name + ".tmp", reference)

new_alignment = idx.get_spliced([0], [genomeSize], strand = "+1")

AlignIO.write(new_alignment, output_prefix + ".fasta", "fasta")

fastaFile = open(output_prefix + ".fasta", 'r')
nexusFile = open(output_prefix + ".nexus", 'w')

alignment = AlignIO.read(output_prefix + ".fasta", 'fasta',
    alphabet=Gapped(IUPAC.ambiguous_dna, '-'))
AlignIO.write(alignment, output_prefix + ".nexus", 'nexus')

fastaFile.close()
nexusFile.close()
