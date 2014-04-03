import os
import datetime
import sys
import numbers
import pprint
import time
from AdamskiClass import AdamskiClass
import tkinter

def main():
    sys.setrecursionlimit(20000)
    list_dna = AdamskiClass.readAndBuildListFromFile('/home/adam/Downloads/rosalind_hamm.txt')
    print(AdamskiClass.computeHammingDistance(list_dna[0],list_dna[1]))

if __name__ == "__main__":
    main()