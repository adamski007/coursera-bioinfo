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

"""
CODE CHALLENGE: Solve the k-Universal Circular String Problem.
Chapter 57.
"""

list_str = []

start_time = datetime.datetime.now()
list_str = genomeClass.generate_all_binary_string(3)
genomeClass.build_Overlap_Graph(list_str)

end_time = datetime.datetime.now()

genomeClass.read_data_and_Build_Graph('listDNA.txt',1)

genomeClass.build_count_in_out_edge_all_nodes()
nodes_in,nodes_out  = genomeClass.check_unbalanced_nodes()
start_node          = genomeClass.build_balanced_graph(nodes_in,nodes_out)

#pprint.pprint(genomeClass.de_bruijn_grapth)
eul_path            = genomeClass.find_eulerian_path(start_node)
str_genome          = genomeClass.string_reconstruction(eul_path)
print str_genome
print len(eul_path)


#print 'Building all binary string took                      : ',end_time-start_time
"""
print ''
start_time = datetime.datetime.now()
genomeClass.build_DeBruijn_Graph_from_listKmers(list_str)
end_time = datetime.datetime.now()
print 'Construction of De Bruijn graph has been builded in  : ',end_time-start_time
"""
"""
genomeClass.build_count_in_out_edge_all_nodes()
list_nodes_unbalanced_in,list_nodes_unbalanced_out = genomeClass.check_unbalanced_nodes()
print 'is the graph un-balanced, list_nodes_in : ',list_nodes_unbalanced_in
print 'is the graph un-balanced, list_nodes_out: ',list_nodes_unbalanced_out
new_edge = genomeClass.build_balanced_graph(list_nodes_unbalanced_in,list_nodes_unbalanced_out)

print 'The new edge inserted and where we have to start from : ',new_edge
"""
"""
start_time = datetime.datetime.now()
eul_path = genomeClass.find_eulerian_path()
print 'Here is the eulerian path : '
print eul_path
end_time   = datetime.datetime.now()
print 'Search of the eulerian path took                     : ',end_time-start_time
"""
#genomeClass.print_eulerian_path(eul_path,0)
#print eul_path
"""
start_time = datetime.datetime.now()
genome = genomeClass.string_reconstruction(eul_path,0)
end_time   = datetime.datetime.now()
print 'String reconstruciton from eulerian path took        : ',end_time-start_time
print genome
print 'Len of overlap graph : ',len(genomeClass.overlap_graph)
"""
