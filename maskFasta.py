#!/usr/bin/env python

import sys
import os
import argparse


def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Mask genes in a fasta file specified \
by a bed file')
    parser.add_argument("fasta", help="fasta file", type=is_file)
    parser.add_argument("bed", help="bedfile", type=is_file)
    parser.add_argument("-mc", help="masking character (default: '-')")
    return parser.parse_args()

def getHeaders(file):
	"""saves fasta header names"""
	f = open("fastaHeaders.txt", 'w')
	for line in open(file):
 		if ">" in line:
   			f.write(line)
	f.close()

def getChromosome(bed):
    """takes chromosome from bed file"""
    for line in open(bed):
        line = line.strip()
        chromosome = line.split('\t')[0]
        print chromosome
    return chromosome

def replaceHeaders(file, chromosome):
	"""changes fasta header names to match bed file"""
	f = open("headerMatch.fasta", 'w')
	for line in open(file):
		if ">" in line:
   			f.write(">" + chromosome + "\n")
		else:
			f.write(line)
	f.close()

def mask(inn, out, bed, mc='\-'):
	"""run maskFastaFromBed"""
	call = 'maskFastaFromBed -fi ' + inn + ' -bed ' + bed +' -fo '+ out +' -mc ' + mc
	os.system(call)

def revert(inn, out, names):
	"""replaces fasta header names with original"""
	fasta= open(inn)
	newnames= open(names)
	newfasta= open(out, 'w')

	for line in fasta:
    		if line.startswith('>'):
        		newname= newnames.readline()
        		newfasta.write(newname)
    		else:
        		newfasta.write(line)

	fasta.close()
	newnames.close()
	newfasta.close()

args = get_args()

## create the name for the masked output file
split = args.fasta.rsplit('/',1)
outfile = "masked_" + split[len(split)-1]

getHeaders(args.fasta)
chromosome = getChromosome(args.bed)
replaceHeaders(args.fasta, chromosome)
if(args.mc):
	mask("headerMatch.fasta","out.fasta", args.bed, args.mc)
else:
	mask("headerMatch.fasta","out.fasta", args.bed)
revert('out.fasta', outfile, 'fastaHeaders.txt')

## remove temporary files
os.system('rm fastaHeaders.txt')
os.system('rm out.fasta')
os.system('rm headerMatch.fasta')
