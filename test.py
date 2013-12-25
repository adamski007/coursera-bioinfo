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
global_alignement = 1

count,matrix_backtrack = genomeClass.lcs(str_v,str_w,matrix_score_lcs,indel_penalty)
#print count
print count[len(str_v),len(str_w)]
#print ''
#print matrix_backtrack
#print ''
global_alignement = 1
print str_v
genomeClass.output_lcs(matrix_backtrack,str_w,len(str_v),len(str_w),global_alignement)
# Below ---
# Print almost the proper string for str_v [ str 1 ] when the second string is the longer
# may miss some char. ?!?!
global_alignement = 0
print '\n'
genomeClass.output_lcs(matrix_backtrack,str_v,len(str_v),len(str_w),global_alignement)
#print '\n'
#print str_w

# NEED TO HANDLE THE CASE , for the other string, both string can contains - character. in the case
# above, we do the handling for str_w and not str_v.