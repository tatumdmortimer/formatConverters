#!/usr/bin/env python

# This script will take an EMBL annotation file and convert to GFF format for
# input into SNPeff. Also needs core alignment with "reference" genome as first
# entry

import sys
import os
import argparse
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(description='Converts EMBL to GFF format')
    parser.add_argument("embl", help="EMBL file name")
    parser.add_argument("core", help="Core genome alignment in fasta format")
    return parser.parse_args()


def convert_embl(embl, fasta):
    f = list(SeqIO.parse(fasta, "fasta"))
    ref = f[0]
    ref.id = "core_genome"
    align_length = len(ref.seq)
    e = SeqIO.read(embl, "embl")
    with open(os.path.basename(embl) + ".gff", "w") as gff_file:
        gff_file.write(
            "##gff-version 3\n##sequence-region core_genome 1 {0}\n".format(align_length))
        for ann in e.features:
            if ann.location.strand == 1:
                strand = "+"
            elif ann.location.strand == -1:
                strand = "-"
            else:
                strand = "."
            gff_file.write(
                "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(
                    "core_genome", "Roary", "CDS", str(ann.location.start + 1),
                    str(ann.location.end + 1), ".", strand, "0",
                    "ID={0};gene={0};locus_tag={0}".format(
                        ann.qualifiers['label'][0])))
        gff_file.write("##FASTA\n")
        SeqIO.write(ref, gff_file, "fasta")


args=get_args()
convert_embl(args.embl, args.core)
