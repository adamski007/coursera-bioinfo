import os
import sys
import numpy
from AdamskiClass import AdamskiClass


infile = open(sys.argv[1],'r')
#infile = open("data.txt",'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)
#listKmers = genomeClass.mostFrequentKmersXMissMatches2(10,2)


#EX : Find most probable kmers in a list of kmers
"""
matrixproba = AdamskiClass.buildMatrixProfile('matrixproba.txt')
print matrixproba
print genomeClass.findMostProbableKmers(5,matrixproba)
"""