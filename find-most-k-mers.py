#!/usr/bin/python

import os;
import sys;

def main_func(input1,length1):
	inputFile = open(input1,'r');
	str 	=	inputFile.readline();
	str = str.replace('\n','');
	len_kmers	= int(length1);
	idx	=	0;
	list_kmers = [];
	most_present_kmers	= 0;
	max_count	= 0;
	while idx <= ( (len(str)) - len_kmers ):
		str_kmers = str[idx:(idx+len_kmers)];
		count = find_x_kmers(str,str_kmers,len_kmers);
		if ( count > max_count ):
			list_kmers = [];
			list_kmers.append(str_kmers);
			max_count = count
		elif ( count == max_count ):
			if not str_kmers in list_kmers:
				list_kmers.append(str_kmers);
				max_count = count;
		idx = idx + 1;
	print(list_kmers);
	inputFile.close();


def find_x_kmers(str,str_kmers,len_kmers):
	idx = 0;
	count = 0;
	while idx <= ( len(str)-len_kmers):
		if ( str[idx:(idx+len_kmers)] == str_kmers ):
			count = count + 1;	
		idx = idx +1;
	return count;

input1=sys.argv[1];
length1 = sys.argv[2];
main_func(input1,length1);
