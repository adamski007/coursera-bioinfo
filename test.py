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

lengthkmers = 3
listkmers = genomeClass.findStringComposition(lengthkmers)
genomeClass.build_Overlap_Graph(listkmers)
pprint.pprint(genomeClass.overlap_graph)
genomeClass.print_edge_overlap_graph()


print '#############################'
print ''
sys.exit(0)
genomeClass.build_DeBruin_Graph(genomeClass.genome,3)
listkmers = genomeClass.findStringComposition(4)
genomeClass.build_Overlap_Graph(listkmers)
genomeClass.print_edge_overlap_graph()
print ''
print '+++++++++++++++++++++++++'
print ''
print genomeClass.de_bruijn_grapth
genomeClass.print_edge_debruijn_graph()