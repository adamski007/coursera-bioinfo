import os
import datetime
import time
import sys
import numpy
import pprint
from AdamskiClass import AdamskiClass


"""
We first need to set our date this way :
our_custom_date = datetime.date(2025,3,28)
timestamp_file  = os.path.getmtime(path)
Then convert this timestamp float in the format date for executing a comparison :
ts_from_file = datetime.date.fromtimestamp(in_float)

Now we can compare both.

"""

our_custom_date = datetime.date(2014,2,28)
file            = '/home/adam/Documents/coursera-courses/bioinformatics-alog1/genome.txt'
timestamp_file  = os.path.getmtime(file)
ts_file         = datetime.date.fromtimestamp(timestamp_file)

print 'Our custom date : ',our_custom_date
print 'Our timestamp from the file : ', timestamp_file
print 'Convertion of this ts to a date is : ', ts_file

if our_custom_date < ts_file:
    print 'Our custom date is older than the file...'
else:
    print 'Our custom date is __NOT__ older than the file...'

sys.exit()


####################################################
path='/home/adam/Documents/coursera-courses/bioinformatics-alog1/genome.txt'
#path='/home/adam/Documents/coursera-courses/bioinformatics-alog1/'
print os.path.basename(path)
path_cwd = os.getcwd()
print path_cwd
print os.stat(path_cwd)
print os.listdir(path_cwd)


# The timestamp of the file is a float.
print 'Type of getmtime : ',type(os.path.getmtime(path))
in_float = os.path.getmtime(path)
in_date = time.ctime(os.path.getmtime(path))
print 'Both type are below '
print in_float
print in_date

# Building a date of choice.
# type of datetime, and not a float
d1 = datetime.date(2025,3,28)
print d1.ctime()
print d1.toordinal()
d1.ctime()
our_datum = d1.ctime()
print 'next statement are printed the same way...'
print datetime.date.fromtimestamp(in_float)
print d1

ts_from_file = datetime.date.fromtimestamp(in_float)
ts_from_our_date  = d1

if ts_from_file > ts_from_our_date:
    print 'youhou working'


if our_datum < in_date:
    print "Our datume is smaller : ",our_datum
else:
    print 'our datume is NOT smaller'


###################################################

sys.exit(0)
sys.setrecursionlimit(20000)
infile = open(sys.argv[1],'r')
#infile = open("data.txt",'r')
nucleotide = infile.readline()
nucleotide = nucleotide.replace('\n','')
genomeClass = AdamskiClass(nucleotide)


list_str = AdamskiClass.readAndBuildListFromFile('listDNA.txt')


"""
list = genomeClass.fill_matrix_BWT(list_str[0])
new_list = ['a','b','a','c','a','d']
print genomeClass.get_x_elem(new_list,'a',3)
"""