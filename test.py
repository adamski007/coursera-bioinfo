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

list_str = AdamskiClass.readAndBuildListPermutation('listDNA.txt')
print list_str[0]
print list_str[1]
print list_str[421]
print len(list_str)

a=[1,2,3]
print genomeClass.sorting_reversal(a)
sys.exit()
print str_v
str_w = list_str[1]
print str_w
str_x = list_str[2]
print str_x