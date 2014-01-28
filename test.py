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
pprint.pprint(bw[0])
last_column = bw[len(bw[0])-1]
last_to_first = genomeClass.last_to_first(bw)
pattern = 'CCT'
print genomeClass.BW_Matching([], last_column, pattern, last_to_first)