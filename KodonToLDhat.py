#!/usr/bin/env python

import sys
import re
import itertools

#determines whether a SNP is valid, only removes N
def snpValidity_N(SNP):
    isValid = False
    no_N = re.compile("N")
    if no_N.match(SNP) is None:
        isValid = True
    return isValid

def snpValidity_all(SNP):
    isValid = False
    no_ambi = re.compile(r"[NYWSKMRDVH]")
    if no_ambi.match(SNP) is None:
        isValid = True
    return isValid

#splits sequences into a defined length
def split_len(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]

#converts MUMmer SNP table to sites.txt and locs.txt for LDhat

#check for correct command line arguments
if len(sys.argv) != 5 :
    print("Usage: KodonToLDhat.py <inputfile> <outputfile prefix> <N or ALL> \
           <name file (ref sequence first)>")
    sys.exit(0)

infilename = sys.argv[1]    #Kodon SNP table
outfilename = sys.argv[2]   #Prefix for sites.txt and locs.txt output files
RefName = sys.argv[3]       #Name of reference sequence
NorALL = sys.argv[4]        #N removed or all ambiguous sequences

InFile = open(infilename, 'r')
outfile_sites = open(outfilename+"_sites.txt", 'w')
outfile_locs = open(outfilename+"_locs.txt", 'w')
PositionList = []
RefSeq = []
NameList = []
SeqList = []
linenumber = 0

#checks validity of SNPS, if valid adds the postion to position list and
#nucleotide to reference sequence

#make a list of all the valid SNPs as lists, [position, SNP, sequence name], in the file
for Line in InFile:
    if linenumber > 1 :	
        Line = Line.strip('\n')
        Line = Line.upper()
        WordList = Line.split('\t')
        SNP = WordList[3][3]
	position = WordList[0]
	name = WordList[2]
        #check validity of the SNP
        if NorALL is 'N':
	    if (snpValidity_N(SNP)) :    
                if position not in PositionList:
                    PositionList.append(position)
                    RefSeq.append(WordList[3][0])
                if name not in NameList:
                    NameList.append(name)   
        else:
            if (snpValidity_all(SNP)) :    
                if position not in PositionList:
                    PositionList.append(position)
                    RefSeq.append(WordList[3][0])
                if name not in NameList:
                    NameList.append(name)
    linenumber = linenumber + 1

#creates a dictionary for each name:sequence
NameList.append(RefName)
for name in NameList:
    SeqList.append(RefSeq[:])
SeqDict = dict(itertools.izip(NameList,SeqList))

#edit sequences in dictionary according to snp table
InFile.seek(0)
linenumber = 0
for line in InFile:
    if linenumber > 1:
        line = line.strip('\n')
        line = line.upper()
        wordlist = line.split()
        SNP = wordlist[3][3]
        seqName = wordlist[2]
        position = wordlist[0]
        if NorALL is 'N':
            if (snpValidity_N(SNP)) and line.find(seqName) != -1:
                pos = PositionList.index(position)
                SeqDict[seqName][pos] = SNP
        else:
            if (snpValidity_all(SNP)) and line.find(seqName) != -1:
                pos = PositionList.index(position)
                SeqDict[seqName][pos] = SNP
    linenumber = linenumber + 1

#write file header for sites.txt
outfile_sites.write(str(len(NameList)) + " " + str(len(PositionList)) + " 1\n")
#write sites.txt
for key in SeqDict:
    outfile_sites.write(">" + key + "\n")
    for i in split_len(''.join(SeqDict[key]), 2000):
        outfile_sites.write(i + "\n")

#write file header for locs.txt
outfile_locs.write(str(len(PositionList)) + " " + PositionList[-1] + " C\n")
#write positions into locs.txt
for position in PositionList:
    outfile_locs.write(position + " ")

InFile.close()
outfile_sites.close()
outfile_locs.close()
