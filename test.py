import os
import datetime
import time
import sys
import numpy
import pprint
from AdamskiClass import AdamskiClass


sys.setrecursionlimit(20000)
infile = open(sys.argv[1],'r')
#infile = open("data.txt",'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)

print genomeClass.findClump_new(5,517,17)

# list_str = AdamskiClass.readAndBuildListFromFile('listDNA.txt')