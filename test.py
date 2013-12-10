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
#eul_path = genomeClass.find_eulerian_path()
#genomeClass.print_eulerian_path(eul_path)
genomeClass.build_count_in_out_edge_all_nodes()

list_nodes_unbalanced_in,list_nodes_unbalanced_out = genomeClass.check_unbalanced_nodes()

pprint.pprint(genomeClass.de_bruijn_grapth)
genomeClass.build_balanced_graph(list_nodes_unbalanced_in,list_nodes_unbalanced_out)
pprint.pprint(genomeClass.de_bruijn_grapth)

#print list_nodes_unbalanced_in
#print list_nodes_unbalanced_out
#print pprint.pprint(genomeClass.de_bruijn_grapth)
#path = genomeClass.find_eulerian_path()
#print path
#genomeClass.print_eulerian_path(path)