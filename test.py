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


listSubStr = genomeClass.findStringComposition(100)
for elem in listSubStr:
    print elem


"""
EX : Find most probable kmers in a list of kmers
matrixproba = AdamskiClass.buildMatrixProfile('matrixproba.txt')
print genomeClass.findMostProbableKmers(6,matrixproba)
"""

"""
FIND MOTIF EXERCICE
listDNA = genomeClass.readAndBuildListFromFile("listDNA.txt")
print AdamskiClass.findMedianString(listDNA,6)
"""

"""
MOTIF ENUMERATION WORKING
listDNA = genomeClass.readAndBuildListFromFile("listDNA.txt")
print listDNA
allmotif = genomeClass.motifEnumeration(listDNA,15,4)
print allmotif.keys()
"""