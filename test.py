import os
import sys
from AdamskiClass import AdamskiClass


infile = open(sys.argv[1],'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)
#complement = genomeClass.buildReverseComplement()
size_kmers = 11
size_window = 566
num_kmers = 18
kmers = "GAAATACTATA"
xMissMatches = 6
list = genomeClass.computeSkewGC()
out = genomeClass.getIdxMinSkewGC(list)
print list
print out
#list = genomeClass.findClump(size_window,size_kmers,num_kmers)


