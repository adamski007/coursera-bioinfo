import os
import datetime
import sys
import numbers
import pprint
import time
from AdamskiClass import AdamskiClass

def main():
    sys.setrecursionlimit(20000)
    # list_dna = AdamskiClass.readAndBuildListFromFile('/home/adam/Downloads/rosalind_hamm.txt')
    # print(AdamskiClass.computeHammingDistance(list_dna[0],list_dna[1]))
    dic_fasta = AdamskiClass.load_fasta_files('test-fasta-file.txt')
    AdamskiClass.get_consensus_profile( dic_fasta)
    """
    list_dna = dic_fasta.values()
    matrix_count = AdamskiClass.get_Count_Matrix_Motifs(list_dna)
    matrix_transposed = matrix_count.transpose()
    print(matrix_transposed.argmax(0))
    print(matrix_transposed[0])
    """

if __name__ == "__main__":
    main()