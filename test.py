import os
import sys
from AdamskiClass import AdamskiClass


infile = open(sys.argv[1],'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)
#complement = genomeClass.buildReverseComplement()
length = 13
list = genomeClass.findMostFrequentKMers(length)
print list
