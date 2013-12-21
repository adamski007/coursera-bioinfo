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
genomeClass.init_start_node()
pprint.pprint(genomeClass.edges_weigth_dag)
print 'Here is the predecessor structure : '
pprint.pprint(genomeClass.predecessor_nodes)
print '_'
#print genomeClass.predecessor_nodes['1']
#print genomeClass.edges_weigth_dag['0'][0][0]
start_node = '10'
another_node = '30'
#genomeClass.weight_nodes[start_node] = 0
test_node = '21'
#print 'Node : ',test_node, 'and his weight is : ', genomeClass.get_weight_nodes(test_node)
print 'Node : ',another_node, 'and his weight is :  ', genomeClass.get_weight_nodes(another_node)
print genomeClass.node_used
print ''
print 'all graph'
pprint.pprint(genomeClass.weight_nodes)
