import os
import datetime
import sys
import numpy
import pprint
import time
from AdamskiClass import AdamskiClass

sys.setrecursionlimit(20000)
infile = open(sys.argv[1],'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)

dic_dna = AdamskiClass.load_fasta_files( sys.argv[1])
gc_max = -1
fasta_id_max_gc = ''
for id_fasta in dic_dna.keys():
    cur_gc = AdamskiClass.compute_GC_content(dic_dna[id_fasta])
    if cur_gc > gc_max:
        fasta_id_max_gc = id_fasta
        gc_max = cur_gc
print fasta_id_max_gc
print gc_max


list_str = AdamskiClass.readAndBuildListFromFile('listDNA.txt')