import os
import sys
from AdamskiClass import AdamskiClass


infile = open(sys.argv[1],'r')
#infile = open("data.txt",'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)
#listKmers = genomeClass.mostFrequentKmersXMissMatches2(10,2)

listDNA = genomeClass.readAndBuildListFromFile("listDNA.txt")
print listDNA
allmotif = genomeClass.motifEnumeration(listDNA,15,4)
print allmotif.keys()