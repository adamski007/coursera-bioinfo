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


"""
list = genomeClass.fill_matrix_BWT(list_str[0])
new_list = ['a','b','a','c','a','d']
print genomeClass.get_x_elem(new_list,'a',3)
"""