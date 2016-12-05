#!/usr/bin/env python

# This script will take a bed file of intervals and an alignment and produce
# a separate fasta alignment roduce nexus alignments of fragments > 1kb

import sys
import argparse
from Bio import AlignIO


def get_args():
    parser = argparse.ArgumentParser(description='Produces alignments from bed \
file intervals and fasta alignment')
    parser.add_argument("align", help="Fasta alignment",
                        type=argparse.FileType('r'))
    parser.add_argument("bed", help="BED file with intervals")
    return parser.parse_args()


def read_bed_file(bed):
    intervalDict = {}
    bedfile = open(bed, "r")
    for index, line in enumerate(bedfile):
        if index > 0:
            line = line.strip().split()
            start = int(line[1])
            stop = int(line[2])
            name = line[3]
            orientation = line[5]
            intervalDict[name] = (start, stop, orientation)
    bedfile.close()
    return intervalDict


def make_alignments(intervalDict, align):
    alignment = AlignIO.read(align, "fasta")
    for interval, startStop in intervalDict.iteritems():
        newAlign = alignment[:, startStop[0]:startStop[1]]
        if startStop[2] == "-":
            for rec in newAlign:
                rec.seq = rec.seq.reverse_complement()
        AlignIO.write(newAlign, interval + ".fasta", 'fasta')

args = get_args()
intervalDict = read_bed_file(args.bed)
make_alignments(intervalDict, args.align)
