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

genomeClass.read_data_and_Build_Graph('listDNA.txt')

print pprint.pprint(genomeClass.de_bruijn_grapth)
print ''
genomeClass.build_count_in_out_edge_all_nodes()
print pprint.pprint(genomeClass.count_in_out_edges)
print ''
list_nodes_unbalanced = genomeClass.check_unbalanced_nodes()
print list_nodes_unbalanced
#print pprint.pprint(genomeClass.de_bruijn_grapth)
#path = genomeClass.find_eulerian_path()
#print path
#genomeClass.print_eulerian_path(path)