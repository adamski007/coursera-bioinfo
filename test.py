import os
import sys
import numpy
from AdamskiClass import AdamskiClass


infile = open(sys.argv[1],'r')
#infile = open("data.txt",'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)

listDNA = AdamskiClass.readAndBuildListFromFile('listDNA.txt')
listMotif = genomeClass.greedy_Motif_Search(listDNA,12,25)
for elem in listMotif:
    print elem