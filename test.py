import os
import datetime
import sys
import numpy
import pprint
import time
from AdamskiClass import AdamskiClass

def main():
    sys.setrecursionlimit(20000)
    infile = open(sys.argv[1],'r')
    nucleotide = infile.readline()
    nucleotide = nucleotide.replace('\n','')
    genomeClass = AdamskiClass(nucleotide)
    genomeClass.loadTableSpectrum(sys.argv[2])
    list = genomeClass.computeMassSpectrum( [genomeClass.genome])
    print list[1]


list_str = AdamskiClass.readAndBuildListFromFile('listDNA.txt')

if __name__ == "__main__":
    main()