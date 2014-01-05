import os
import datetime
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

list_str = AdamskiClass.readAndBuildListFromFile('listDNA.txt')
str_v = list_str[0]
print str_v
str_w = list_str[1]
print str_w
str_x = list_str[2]
print str_x


genomeClass.score_multiple_alignement(list_str)
print genomeClass.matrix_multiple_alignment.shape
print genomeClass.matrix_multiple_alignment
print genomeClass.matrix_multiple_alignment[len(str_v)][len(str_w)][len(str_x)]

sys.exit()
#print count
print int(count[len(str_v),len(str_w)])
global_alignement = True
first_str = True
genomeClass.output_lcs(matrix_backtrack,str_v,len(str_v),len(str_w),global_alignement,first_str)
print ''
first_str = False
genomeClass.output_lcs(matrix_backtrack,str_w,len(str_v),len(str_w),global_alignement,first_str)