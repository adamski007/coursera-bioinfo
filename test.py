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
str_v = list_str[0]
str_w = list_str[1]

matrix_score_lcs = genomeClass.build_scoring_matrix('score_matrix.txt')
indel_penalty = 5

local_alignment=True
count,matrix_backtrack = genomeClass.lcs(str_v,str_w,matrix_score_lcs,indel_penalty,local_alignment)
print count.max()
#print matrix_backtrack
print int(count[len(str_v),len(str_w)])
global_alignement = True
first_str = True
genomeClass.output_lcs(matrix_backtrack,str_v,len(str_v),len(str_w),global_alignement,first_str)
print ''
first_str = False
genomeClass.output_lcs(matrix_backtrack,str_w,len(str_v),len(str_w),global_alignement,first_str)

print ''
print count.shape

idx_i = genomeClass.max_value_matrix_score[0]
idx_j = genomeClass.max_value_matrix_score[1]
print 'Idx i and j ',idx_i,idx_j
print genomeClass.backtrack_local_alignement(count,idx_i,idx_j,str_v,str_w,indel_penalty,matrix_score_lcs)