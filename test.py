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
amino = 'NQEL'
listSubPep = genomeClass.generateListSubPep(amino)
listMassSpectrum = genomeClass.computeMassSpectrum(listSubPep)
spectrum="1189 137 771 245 833 234 71 842 1035 454 186 461 202 980 330 503 216 1172 970 743 0 640 663 1069 532 1058 584 1087 358 447 113 1190 794 711 470 719 131 1101 525 268 1200 986 461 874 956 663 931 420 849 202 1172 131 114 973 1101 372 697 103 429 323 1232 770 987 592 800 1117 1166 945 606 1117 656 560 347 509 1303 186 1076 317 533 333 316 339 856 883 647 227 131 778 964 842 1172 64"
spectrumList    =   spectrum.split(' ')
spectrumListInt = []
for elem in spectrumList:
    spectrumListInt.append(int(elem))
spectrumListInt.sort()
listResult = genomeClass.computeSpectralConvolution(spectrumListInt)
listResult.sort()
print listResult
