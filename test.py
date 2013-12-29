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

matrix_score_lcs = genomeClass.build_scoring_matrix('score_matrix_fit_align.txt')
indel_penalty = int(1)

edit_distance = False
local_alignment = False
fitting_alignment = True
count,matrix_backtrack = genomeClass.lcs(str_v,str_w,matrix_score_lcs,indel_penalty,local_alignment,edit_distance,fitting_alignment)
idx_i , idx_j = AdamskiClass.get_idx_max_value_last_column_matrix(count)
print int(count[idx_i,idx_j])


global_alignement = True
first_str = True
fitting_alignment = True
genomeClass.output_lcs(matrix_backtrack,str_v,idx_i,idx_j,global_alignement,first_str,fitting_alignment)
print ''
first_str = False
genomeClass.output_lcs(matrix_backtrack,str_w,idx_i,idx_j,global_alignement,first_str,fitting_alignment)