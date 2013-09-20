#!/usr/bin/env python

import sys
import re
import itertools

  
#SNP table from VCFstoSNPTable.py to LDhat files

#determines whether a SNP is valid, only removes N
def snpValidity_N(SNP) :
    isValid = False
    no_N = re.compile("N")
    if no_N.match(SNP) is None:
        isValid = True
    return isValid

#determines whether a SNP is valid, removes all ambiguous bases
def snpValidity_all(SNP) :
    isValid = False
    no_ambi = re.compile(r"[NYWSKMRDVH]")
    if no_ambi.match(SNP) is None:
        isValid = True
    return isValid

#splits sequences into a defined length
def split_len(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]



#check for correct commandline arguments
if len(sys.argv) != 5 :
	print("Usage:  SNPTableToLDhat.py  <input file>  <outputprefix> \
               <name of reference sequence> <N or ALL>")
	sys.exit(0)

InFileName = sys.argv[1]   #SNP table
OutFileName = sys.argv[2]  #Nexus format
Refname = sys.argv[3] #Reference sequence
NorALL = sys.argv[4] #Should all ambiguous bases be removed or just N

InFile = open(InFileName, 'r')	#access the file to be read from
outfile_sites = open(OutFileName + '_sites.txt', 'w') 
outfile_locs = open(OutFileName + '_locs.txt', 'w') 
PositionList = [] 
NameList = []
Refseq = []
SeqList = []
linenumber = 0

for Line in InFile:
    Line = Line.strip('\n')
    WordList = Line.split()
    SNP = WordList[3][3]
    position = WordList[2]
    name = WordList[0]
    #check validity of the SNP
    if NorALL is 'N':
        if (snpValidity_N(SNP)) :    
            if position not in PositionList:
                PositionList.append(position)
                Refseq.append(WordList[3][0])
            if name not in NameList:
                NameList.append(name)   
    else:
        if (snpValidity_all(SNP)) :    
            if position not in PositionList:
                PositionList.append(position)
                Refseq.append(WordList[3][0])
            if name not in NameList:
                NameList.append(name)
print PositionList
#creates a dictionary for each name:sequence
for name in NameList:
    SeqList.append(Refseq[:])
SeqDict = dict(itertools.izip(NameList,SeqList))

#edit sequences in dictionary according to snp table
InFile.seek(0)
for line in InFile:
    line = line.strip('\n')
    wordlist = line.split()
    SNP = wordlist[3][3]
    seqName = wordlist[0]
    position = wordlist[2]
    if NorALL is 'N':
        if (snpValidity_N(SNP)) and line.find(seqName) != -1:
            pos = PositionList.index(position)
            SeqDict[seqName][pos] = SNP
    else:
        if (snpValidity_all(SNP)) and line.find(seqName) != -1:
            pos = PositionList.index(position)
            SeqDict[seqName][pos] = SNP


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
