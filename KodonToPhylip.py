#!/usr/bin/env python

import sys
import re

#converts Kodon alignment into a phylip file

#searches file for lines that begin with a particular sequence name
def write_sequence(seq_name, inputfile):
    inputfile.seek(0)
    sequence = ""    
    for line in inputfile:
        line = line.strip('\n')
        if seq_name in line:
            sequence = sequence + line
    sequence = sequence.replace(seq_name,"")
    sequence = sequence.replace(" ","")
    sequence = re.sub(r'[.0-9]', '', sequence)
    sequence = sequence.upper()
    return sequence
             
#check for correct command line arguments
if len(sys.argv) != 6:
    print("Usage: KodonToPhylip.py <inputfile> <outputfile> <namefile> \
           <number of sequences> <number of nucleotides in alignment>")
    sys.exit(0)

infilename = sys.argv[1]    #Kodon alignment file
outfilename = sys.argv[2]   #Phylip format
seq_names = sys.argv[3]     #File with names of sequence in Kodon
seq_num = sys.argv[4]       #number of sequences in alignment
nuc_num = sys.argv[5]       #number of nucleotides in aligned sequences


infile = open(infilename, 'r')
outfile = open(outfilename, 'w')
namefile = open(seq_names, 'r')

#write file header
outfile.write(seq_num + " " + nuc_num + "\n")

#write sequences into file
for line in namefile:
    name_kodon = line.strip('\n')
    name_phylip = name_kodon[:]
    if len(name_phylip) < 10:
        name_phylip = name_phylip + '          '
    name_phylip = name_phylip[:10]
    outfile.write(name_phylip+write_sequence(name_kodon, infile)+'\n')

infile.close()
outfile.close()
namefile.close()
