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
#print 'suite'
genomeClass.build_count_in_out_edge_all_nodes()
list_nodes_unbalanced_in,list_nodes_unbalanced_out = genomeClass.check_unbalanced_nodes()
new_edge = genomeClass.build_balanced_graph(list_nodes_unbalanced_in,list_nodes_unbalanced_out)

print 'The new edge inserted and where we have to start from : ',new_edge
eul_path = genomeClass.find_eulerian_path(new_edge)
genomeClass.print_eulerian_path(eul_path)
print eul_path
genome = genomeClass.string_reconstruction(eul_path)
print genome
bin_str = genomeClass.generate_all_binary_string(4)
print bin_str
new_str = '11'
print genomeClass.add_necessary_zeros(new_str,4)
print new_str
list_str = genomeClass.generate_all_binary_string(18)
print list_str
#genomeClass.print_edge_debruijn_graph()
"""
new_edge = genomeClass.build_balanced_graph(list_nodes_unbalanced_in,list_nodes_unbalanced_out)
print 'Starting node found in building balanced graph : ',new_edge
eul_path = genomeClass.find_eulerian_path(new_edge)
print 'Starting node is : ',new_edge
genomeClass.print_eulerian_path(eul_path)
"""

#print list_nodes_unbalanced_in
#print list_nodes_unbalanced_out
#print pprint.pprint(genomeClass.de_bruijn_grapth)
#path = genomeClass.find_eulerian_path()
#print path
#genomeClass.print_eulerian_path(path)