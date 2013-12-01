import os
import sys
from AdamskiClass import AdamskiClass


infile = open(sys.argv[1],'r')
#infile = open("data.txt",'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)
#listKmers = genomeClass.mostFrequentKmersXMissMatches2(10,2)
#print listKmers
#kmres = 'AAAA'
#listBuilded = []
#listKmers = genomeClass.generateAllKmersXMissMatches(kmres,2,listBuilded)
#print listKmers

x='CCAGGTATCACATAGACACCGCACAAGGATACAGAGAAGCGC'
y='AACCAGAAACGCGTAAATCCACGGTAGCGTGAGTGATAACGC'
w='AAGCGCACTTTAACCTGGGAGGCTTGGGTGCTAGCCCAGAAC'
z='ACCTGAAAACGCCGCACACTATGAGGGTACGCGAGGTCTGTC'
v='TAGTTTGTCCGTAGATGTAAGCGCAGGCCAAATTTCCTGTTT'
listDNA = []
listDNA.append(x)
listDNA.append(y)
listDNA.append(w)
listDNA.append(z)
listDNA.append(v)

kmers='AAG'
DNA='GCAATCCTCAGC'
diff,motif =  AdamskiClass.findMotif(kmers,DNA)
print 'The hammind distance is : ',diff, 'and the kers with this count is : ',motif

total,dna = AdamskiClass.computeSumHammingDistance('AAA',listDNA)
print total,dna

print AdamskiClass.findMedianString(listDNA,3)
listDNA2 = AdamskiClass.readAndBuildListFromFile('listDNA.txt')
print AdamskiClass.findMedianString(listDNA2,6)
