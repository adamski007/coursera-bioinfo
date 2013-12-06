import os
import sys
import numpy
import pprint
from AdamskiClass import AdamskiClass


infile = open(sys.argv[1],'r')
#infile = open("data.txt",'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)

list_kmers = AdamskiClass.readAndBuildListFromFile('listDNA.txt')
#genomeClass.build_Overlap_Graph(list_kmers)
#genomeClass.print_edge_overlap_graph()
#print ''
genomeClass.build_DeBruijn_Graph_from_listKmers(list_kmers)
genomeClass.print_edge_debruijn_graph()