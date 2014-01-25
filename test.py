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
#print 'Starting programme : ',time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
print genomeClass.get_inverse_BW_transform(list_str[0])
#print genomeClass.get_x_elem(new_list,'a',3)