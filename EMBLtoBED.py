#!/usr/bin/env python

# This script will take an EMBL annotation file and convert to BED format

import sys
import os
import argparse
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(description='Converts EMBL to BED format')
    parser.add_argument("embl", help="EMBL file name")
    return parser.parse_args()


def convert_embl(embl):
    e = SeqIO.read(embl, "embl")
    with open(os.path.basename(embl) + ".bed", "w") as bedfile:
        for ann in e.features:
            if ann.location.strand == 1:
                strand = "+"
            elif ann.location.strand == -1:
                strand = "-"
            else:
                strand = "."
            bedfile.write(
                "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                    e.id, str(ann.location.start),
                    str(ann.location.end),
                    ann.qualifiers['label'][0],
                    '1000', strand))


args = get_args()
convert_embl(args.embl)
