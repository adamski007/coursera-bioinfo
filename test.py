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

genomeClass.read_data_weigth_edges('listDNA.txt')
pprint.pprint(genomeClass.edges_weigth_dag)
print genomeClass.edges_weigth_dag['0'][0][0]
