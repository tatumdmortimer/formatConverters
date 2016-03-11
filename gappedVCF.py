#!/usr/bin/env python

import sys
import os
import argparse
import random
from Bio import AlignIO

# Removes sites with gaps from a VCF produced by snp-sites using the alignment

def get_args():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Subsamples a VCF to create \
SFS')
    parser.add_argument("vcf", help="VCF to be sampled",
                    type=argparse.FileType('r'))
    parser.add_argument("align", help="Alignment to find gaps",
                    type=argparse.FileType('r'))
    return parser.parse_args()
    
def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

def edit_vcf(vcf, gapD):
    outVCF = open("gapped.vcf", "w")
    for line in vcf:
        if line[0] == "#":
            outVCF.write(line)
            if line[1:6] == "CHROM":
                d = {}
                headers = line.strip().split()
                for i, samp in enumerate(headers):
                    d[samp] = i
        else:
            line_items = line.strip().split()
            pos = int(line_items[1]) - 1
            if pos in gapD:
                print(str(pos) + '\t' + ",".join(gapD[pos]))
                for s in gapD[pos]:
                    ind = d[s]
                    print("Changing {0} to a gap".format(ind))
                    print("\t".join(line_items))
                    line_items[ind] = "-"
                    print("\t".join(line_items))
            outVCF.write("\t".join(line_items) + '\n')
    vcf.close()
    outVCF.close()

def find_gaps(align):
    alignment = AlignIO.read(align, "fasta")
    gapD = {}
    for seqRecord in alignment:
            gapIndex = find(seqRecord, '-')
            for gap in gapIndex:
                if gap in gapD:
                    gapD[gap].append(seqRecord.id)
                else:
                    gapD[gap] = [seqRecord.id]
    return gapD

args = get_args()

gapD = find_gaps(args.align)
edit_vcf(args.vcf, gapD)
