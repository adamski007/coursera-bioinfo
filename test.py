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
genomeClass.loadTableSpectrum(sys.argv[3])
#print genomeClass.codonTable['UGC']
amino = 'YMRPTMNVQEPEAFP'
listSubPep = genomeClass.generateListSubPep(amino)
listMassSpectrum = genomeClass.computeMassSpectrum(listSubPep)
print listMassSpectrum