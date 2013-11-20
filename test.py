import os
import sys
from AdamskiClass import AdamskiClass


infile = open(sys.argv[1],'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)
#complement = genomeClass.buildReverseComplement()
size_kmers = 11
size_window = 566
num_kmers = 18
kmers = "GAAATACTATA"
xMissMatches = 6
genomeClass.loadCodonTable(sys.argv[2])
#print genomeClass.codonTable['UGC']
acido = 'ATTACCAAGGCAAATATCCCCATAAAACTATTTCCC'
x = genomeClass.findSubStringAcido(acido)
print x,


"""
print genomeClass.genome
rna = genomeClass.transcribeDNAToRNA()
print genomeClass.translateRNAIntoAcido(rna)
reverseComplement = genomeClass.buildReverseComplement()
print reverseComplement
genomeClass.genome = reverseComplement
rna = genomeClass.transcribeDNAToRNA()
print genomeClass.translateRNAIntoAcido(rna)
"""