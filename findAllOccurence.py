#!/usr/bin/python

import os;
import sys;

def findAllOcc(inputStr,kmers):
	listPosition = []
	idx = 0
	while idx <= ( len(inputStr) - len(kmers) ):
		if ( kmers == inputStr[idx:idx+len(kmers)] ) :
			listPosition.append(idx)
		idx = idx + 1
	return listPosition

inFile = open(sys.argv[1],'r')
kmers	= sys.argv[2];
genome = inFile.readline()
genome = genome.replace('\n','')
print(findAllOcc(genome,kmers))
