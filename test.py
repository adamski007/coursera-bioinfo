import os
import datetime
import sys
import numpy
import pprint
from AdamskiClass import AdamskiClass


infile = open(sys.argv[1],'r')
#infile = open("data.txt",'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)

matrix_down, matrix_right = AdamskiClass.initialize_matrix_edges('matrix.txt',11,10)
print matrix_down
print ''
print matrix_right
print matrix_right[2][3]
print AdamskiClass.manathan_tourist(11,10,matrix_down,matrix_right)