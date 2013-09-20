#!/usr/bin/env python

import sys
import re

#This script converts a nexus file to phylip file

#check for correct command line arguments
if len(sys.argv) != 3:
    print("Usage: NexusToPhylip.py <inputfile> <outputfile>")
    sys.exit(0)

infilename = sys.argv[1]    #Nexus input file
outfilename = sys.argv[2]   #Phylip format

infile = open(infilename, 'r')
outfile = open(outfilename, 'w')

linenumber = 0
seq_num = ''
nuc_num = ''
for line in infile:
    #write file header
    if linenumber == 3:
        line = line.strip('\n')
        wordlist = line.split()
        seq_num = re.sub("\D", "", wordlist[1])
        nuc_num = re.sub("\D", "", wordlist[2])
        outfile.write(seq_num + " " + nuc_num + ('\n'))
    #write sequences to file
    if linenumber > 5:
        line = line.strip('\n')
        if linenumber % 2 == 0:
            name = line
            if len(name) > 9:
                outfile.write(name[:11])
            else:
                name = name + '          '
                outfile.write(name[:11])
        else:
            outfile.write(line + '\n')
    linenumber = linenumber + 1

infile.close()
outfile.close()
