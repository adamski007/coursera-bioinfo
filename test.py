import os
import datetime
import sys
import numpy
import pprint
from AdamskiClass import AdamskiClass


sys.setrecursionlimit(20000)
infile = open(sys.argv[1],'r')
#infile = open("data.txt",'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)


list_str = AdamskiClass.readAndBuildListFromFile('listDNA.txt')
# list_str will contain only one string in this exercice.
text = list_str[0]

print 'Beginning getting all suffixes @ : ',str(datetime.datetime.now())
list_suffixes = genomeClass.get_all_suffix(text)
print 'All suffixes builded @ : ',str(datetime.datetime.now())
print ''
print '     Beginning construction of suffixes trie structure @ : ',str(datetime.datetime.now())
genomeClass.build_tries_construction(list_suffixes)
print '     Construction of suffixes trie structure finished  @ : ',str(datetime.datetime.now())

print '         Searching for the longest repeat in the trie strucutre @ : ',str(datetime.datetime.now())
genomeClass.search_longest_repeat()
print '         Longest repeat found in the trie strucutre @ : ',str(datetime.datetime.now())
print genomeClass.longest_repeat
