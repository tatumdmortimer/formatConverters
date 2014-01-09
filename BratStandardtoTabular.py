#!/usr/bin/env python

import sys

#This script takes a BratNextGen standard output file and converts it to the
#tabular form

# check for correct arguments
if len(sys.argv) != 3:
    print("Usage: BratStandardtoTabular.py <inputfile> <outputfile>")
    sys.exit(0)

input_name = sys.argv[1]
output_name = sys.argv[2]

input_file = open(input_name, 'r')
output_file = open(output_name, 'w')

output_file.write('Start\tEnd\tOrigin\tStrain\n')
strain = ''
for index,line in enumerate(input_file):
    if index > 1:
        line = line.strip()
        if line[-1] == ':':
            wordList = line.split()
            strain = wordList[2]
        else:
            wordList = line.split()
            output_file.write(wordList[0] + '\t' + wordList[1] + '\t' + wordList[2] + '\t' + strain + '\n')
input_file.close()
output_file.close()
