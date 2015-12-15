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

def edit_vcf(vcf, gapSet):
    outVCF = open("gapless.vcf", "w")
    for line in vcf:
        if line[0] == "#":
            outVCF.write(line)
            continue
        line_items = line.strip().split()
        pos = int(line_items[1]) - 1
        if pos not in gapSet:
            outVCF.write(line)
    vcf.close()
    outVCF.close()

def find_gaps(align):
    alignment = AlignIO.read(align, "fasta")
    gapSet = set()
    for i,seqRecord in enumerate(alignment):
        gapIndex = find(seqRecord.seq, '-')
        for gap in gapIndex:
            gapSet.add(gap)
    return gapSet

args = get_args()

gapSet = find_gaps(args.align)
edit_vcf(args.vcf, gapSet)
