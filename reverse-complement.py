#!/usr/bin/python
import os;
import sys;

def give_complement(letter):
	if ( letter == 'A' ):
		return 'T';
	elif ( letter == 'T' ):
		return 'A';
	elif ( letter == 'C' ):
		return 'G'
	elif ( letter == 'G' ):
		return 'C'

def buildReverseComplement(input):
	complement = ''
	for char in input:
		complement = complement + give_complement(char)
	return complement[::-1]

infile = open(sys.argv[1],'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
print(buildReverseComplement(nucleotide))
infile.close()
