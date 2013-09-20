#!/usr/bin/env python

import sys
import re
import itertools

  
#SNP table from VCFstoSNPTable.py to nexus file 

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



#check for correct commandline arguments
if len(sys.argv) != 5 :
	print("Usage:  SNPTableToNexus.py  <input file>  <outputfile> \
               <name of reference sequence> <N or ALL>")
	sys.exit(0)

InFileName = sys.argv[1]   #SNP table
OutFileName = sys.argv[2]  #Nexus format
Refname = sys.argv[3] #Reference sequence
NorALL = sys.argv[4] #Should all ambiguous bases be removed or just N

InFile = open(InFileName, 'r')	#access the file to be read from
OutFile = open(OutFileName, 'w') #overwrite what is present in the file
PositionList = [] 
NameList = []
Refseq = []
SeqList = []
linenumber = 0

#make a list of all the valid SNPs as lists, [position, SNP, sequence name], in the file
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

#write file header
OutFile.write("#NEXUS" + "\n\n" + "Begin DATA;" + "\n\t" + "Dimensions ntax=" +
	str(len(NameList)+1) + " nchar=" + str(len(PositionList)) + ";" + "\n\t" + "Format datatype=DNA gap=-;" 
	+ "\n\t" + "Matrix" + "\n")

#write ref sequence into file
OutFile.write(Refname + '\n')
OutFile.write("".join(Refseq))
OutFile.write('\n')

#write all of the SNPs into the new file
for key in SeqDict:
    OutFile.write(key + "\n" + ''.join(SeqDict[key]) + "\n")
OutFile.write(";" + "\n" + "END;")

InFile.close()
OutFile.close()
