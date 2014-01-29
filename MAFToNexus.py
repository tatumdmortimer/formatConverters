#!/usr/bin/env python

#import required libraries
import sys
sys.path.insert(1, "/home/peplab/src/alignio-maf")
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

# get alignments with all sequences and reverse complement those that
# are on incorrect strand
for a in multiple_alignment:
    strainIDs = []
    for seqRecord in a:
        strainIDs.append(seqRecord.id)
    if len(strainIDs) == numSeq:
        refInd = strainIDs.index(reference)
        if a[refInd].annotations["strand"] == "-1":
            newAlignment = []
            for seqRecord in a:
                seq = seqRecord.reverse_complement(id=True, name=True,
                    features=True, annotations=True)
                if seq.annotations["strand"] == "-1":
                    seq.annotations["strand"] = "+1"
                    seq.annotations["start"] = int(seq.annotations["start"]) - int(seq.annotations["size"]) + 1
                elif seq.annotations["strand"] == "+1":
                    seq.annotations["strand"] = "-1"
                    seq.annotations["start"] = int(seq.annotations["start"]) + int(seq.annotations["size"]) - 1
                newAlignment.append(seq)
            refList.append(MultipleSeqAlignment(newAlignment))
        else:
            refList.append(a)

# create new temporary MAF file
AlignIO.write(refList, input_name + ".tmp", "maf")

# make MAF index file
idx = MafIO.MafIndex(reference + ".mafindex", input_name + ".tmp", reference)

new_alignment = idx.get_spliced([0], [genomeSize], strand = "+1")

AlignIO.write(new_alignment, output_name, "fasta")

# take blocks that contain reference sequence and add them to a list
#alignment = []
#for a in multiple_alignment:
#    if a[0].id.startswith(ref):
#        alignment.append((int(a[0].annotations["start"]),a))

# sort the block list by start of ref sequence
#sorted_alignment = sorted(alignment, key=itemgetter(0))

# make a dictionary for genome sequences
#genomeDict = {}

# read blocks and edit genome dictionary
#for block in sorted_alignment:
#    block = block[1] # get rid of tuple format    
#    # get start and stop sites for reference
#    blockStart = int(a[0].annotations["start"])
#    blockStop = blockStart + int(a[0].annotations["size"])
#    for seqRecord in block:
#        strainName = seqRecord.id.split('.')[0]
#        if strainName not in genomeDict:
#            # add sequence to genome dictionary
#            # using '-' as placeholder until actual sequence is added
#            genomeDict[strainName] = ['-']*genomeSize
#        genome = genomeDict[strainName]
#        genome[blockStart:blockStop] = list(seqRecord.seq)
#        genomeDict[strainName] = genome
#        
#outfile.close()
