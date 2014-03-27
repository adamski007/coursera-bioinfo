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
dna_string = infile.readline()
dna_string = dna_string.replace('\n','')
dna = AdamskiClass(dna_string)
print dna.buildReverseComplement()
print 2
print 4
print 5

sys.exit()

list_str = AdamskiClass.readAndBuildListFromFile('listDNA.txt')
partial_suffixe_array = genomeClass.build_suffixe_array(list_str[0], True, 7)
genomeClass.print_idx_suffixe_array( partial_suffixe_array, True)