import os
import sys
from AdamskiClass import AdamskiClass


infile = open(sys.argv[1],'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)
#listKmers = genomeClass.mostFrequentKmersXMissMatches2(10,2)
#print listKmers
kmres = 'AAA'
listKmers = []
listKmers2 = genomeClass.generateAllKmersXMissMatches(kmres,listKmers,2)
#print listKmers

print '++++++++'
print ''
listKmers.sort()
print 'sorted list...'
#print listKmers
print 'Size of the list is : ',len(listKmers)

# GCT GGT TGGT
# GCG CAC ACAC





# This comment should leave after git command...








"""
#complement = genomeClass.buildReverseComplement()
size_kmers = 11
size_window = 566
num_kmers = 18
kmers = "GAAATACTATA"
xMissMatches = 6
genomeClass.loadCodonTable(sys.argv[2])
genomeClass.loadTableSpectrum(sys.argv[3])
#print genomeClass.codonTable['UGC']
#print genomeClass.massTable['T']
if 101 in genomeClass.massTable:
    print "101 present in masstable of the class"
amino = 'NQEL'
listSubPep = genomeClass.generateListSubPep(amino)
listMassSpectrum = genomeClass.computeMassSpectrum(listSubPep)
spectrum="1189 137 771 245 833 234 71 842 1035 454 186 461 202 980 330 503 216 1172 970 743 0 640 663 1069 532 1058 584 1087 358 447 113 1190 794 711 470 719 131 1101 525 268 1200 986 461 874 956 663 931 420 849 202 1172 131 114 973 1101 372 697 103 429 323 1232 770 987 592 800 1117 1166 945 606 1117 656 560 347 509 1303 186 1076 317 533 333 316 339 856 883 647 227 131 778 964 842 1172 64"
spectrum="129 510 464 285 902 593 682 852 639 567 373 131 795 170 163 333 551 886 319 156 129 137 620 884 618 958 260 696 730 301 886 585 113 526 1015 395 430 220 266 714 859 755 878 489 376 749 244 397 448 642 505 0 845 749 57 422 266 771"
spectrum="465 473 998 257 0 385 664 707 147 929 87 450 748 938 998 768 234 722 851 113 700 957 265 284 250 137 317 801 128 820 321 612 956 434 534 621 651 129 421 337 216 699 347 101 464 601 87 563 738 635 386 972 620 851 948 200 156 571 551 522 828 984 514 378 363 484 855 869 835 234 1085 764 230 885"
spectrum="743 934 101 588 156 57 557 261 490 1145 376 103 587 115 274 372 289 989 372 97 646 1042 1032 1042 871 500 211 114 158 769 371 216 645 1044 375 271 716 487 429 655 186 1031 959 660 402 929 876 542 658 770 1030 831 928 748 217 103 269 613 532 761 759 987 500 985 384 773 1048 0 645 283 485 884 558 160 773 862 386 856 314 874 472 113 397 216 499 774 929 673 872 603 1088 273"
spectrum = "0 137 186 323"
genomeClass.buildAllMassValue()
spectrumList    =   spectrum.split(' ')
spectrumListInt = []
for elem in spectrumList:
    spectrumListInt.append(int(elem))
spectrumListInt.sort()
test = genomeClass.findOneMers(spectrumListInt)
"""