import os
import datetime
import sys
import numpy
import pprint
import time
from AdamskiClass import AdamskiClass


sys.setrecursionlimit(20000)
infile = open(sys.argv[1],'r')
#infile = open("data.txt",'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)

list_str = AdamskiClass.readAndBuildListFromFile('listDNA.txt')
bw = genomeClass.fill_matrix_BWT(list_str[0], False, True)
pprint.pprint(bw)
print genomeClass.last_to_first(bw)