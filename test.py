import os
import datetime
import sys
import numpy
import pprint
from AdamskiClass import AdamskiClass


sys.setrecursionlimit(2000)
infile = open(sys.argv[1],'r')
#infile = open("data.txt",'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)

list_str = AdamskiClass.readAndBuildListFromFile('listDNA.txt')
#print list_str[0]
#print ''
#print list_str[1]

str_v = list_str[0]
str_w = list_str[1]

matrix_score_lcs = genomeClass.build_scoring_matrix('score_matrix.txt')
indel_penalty = 5

count,matrix_backtrack = genomeClass.lcs(str_v,str_w,matrix_score_lcs,indel_penalty)
print count
#print ''
print matrix_backtrack
#print ''
AdamskiClass.output_lcs(matrix_backtrack,str_v,len(str_v),len(str_w))