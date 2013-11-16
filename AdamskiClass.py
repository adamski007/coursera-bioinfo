class AdamskiClass:
    """ This is a test for the documentation """
    genome = ""
    codonTable = {}

    def __init__(self, strGenome):
        self.genome = strGenome
        self.codonTable = {}

    """
        Input : 1. a letter in the genome.
        Output: The complement of this letter in the genome.
    """
    def give_complement(self, letter):
        if ( letter == 'A' ):
            return 'T'
        elif ( letter == 'T' ):
            return 'A'
        elif ( letter == 'C' ):
            return 'G'
        elif ( letter == 'G' ):
            return 'C'

    """
        Input : 1. a genome
        Ouput : the reverse complement of this genome.
    """
    def buildReverseComplement(self):
        complement = ''
        for char in self.genome:
            complement = complement + self.give_complement(char)
        return complement[::-1]

    """
        Input : 1. a kmers
        Output: a list with the position where the kmers a located in the genome.
    """
    def findAllOcc(self,kmers):
        listPosition = []
        idx = 0
        while idx <= ( len(self.genome) - len(kmers) ):
            if ( kmers == self.genome[idx:idx+len(kmers)] ) :
                listPosition.append(idx)
            idx = idx + 1
        return listPosition

    """ Input :     1. a genome
                    2. a kmers
        Output:     the numbers of kmers we found in the genome.
    """
    def find_x_kmers(self,str_kmers,len_kmers):
        idx = 0;
        count = 0;
        while idx <= ( len(self.genome)-len_kmers):
            if ( self.genome[idx:(idx+len_kmers)] == str_kmers ):
                count = count + 1;
            idx = idx +1;
        return count;

    """
        Input   :   1. a genome
                    2. length of a kmers
        Output  :   a list containing all the kmers with a specific size in the genome.
    """
    def findMostFrequentKMers(self,length1):
        # For historic reasons...
        len_kmers = length1
        idx	=	0
        list_kmers = []
        most_present_kmers	= 0
        max_count	= 0
        while idx <= ( (len(self.genome)) - len_kmers ):
            str_kmers = self.genome[idx:(idx+len_kmers)]
            count = self.find_x_kmers(str_kmers,len_kmers)
            if ( count > max_count ):
                # Re-init the list as we got more a bigger k-mers in this attempt.
                list_kmers = []
                list_kmers.append(str_kmers)
                max_count = count
            elif ( count == max_count ):
                if not str_kmers in list_kmers:
                    list_kmers.append(str_kmers)
                    # Next line probably not needed. TO DO [ to check actually. ]
                    max_count = count
            idx = idx + 1
        return list_kmers

    """
        Input   :   1. a genome
                    2. a size of window where to find the number of kmers
                    3. the length of the kmers
                    4. the number of kmers we need to find the specified window.
        Output  :   The list of kmers found in the window.
    """
    def findClump(self,size_window,size_kmers,number):
        idx = 0
        list_kmers = []
        globalGenome = self.genome
        while idx <= ( len(globalGenome) - size_window ):
            window = self.genome[idx:(idx+size_window)]
            self.genome = window
            # Doing the computation temporarily
            list_kmers_in_window = self.findMostFrequentKMers(size_kmers)
            # Checking if these kmers are present enough in the window.
            for kmers in list_kmers_in_window:
                if self.find_x_kmers(kmers,len(kmers)) >= number:
                    # It does means that we found at least the requested number of kmers in
                    # the specified window. We got one candidate for the the kmers.
                    if not kmers in list_kmers:
                        list_kmers.append(kmers)
            idx = idx + 1
            self.genome = globalGenome
        return list_kmers

    """
        Input   : 1. a genome
        Output  : a list containing at each idx of the genome, the current skew
            The skew is defined as the defirence between G and C at each indices of the genome.
    """
    def computeSkewGC(self):
        skewGC = [0]
        curSkew = 0
        for char in self.genome:
            if char == 'C':
                curSkew = skewGC[-1]-1
                skewGC.append(curSkew)
            elif char == 'G':
                skewGC.append(skewGC[-1]+1)
            else:
                skewGC.append(skewGC[-1])
        return skewGC

    """
        Input   : the skew of a genome
        Output  : a list containing the index where the skew is the minimum.
    """

    def getIdxMinSkewGC(self,listSkew):
        idxMinSkewGC = []
        # Making a copy because the sort is done in place.
        skewCopy = []
        skewCopy.extend(listSkew)
        skewCopy.sort()
        smalestItem = skewCopy[0]
        count = skewCopy.count(smalestItem)
        idx = 0
        # Re resseting the skewCopy for the next while
        skewCopy = []
        skewCopy.extend(listSkew)
        while idx < count:
            # Getting all the indices where the skew is the smallest.
            idxItem = skewCopy.index(smalestItem)
            idxMinSkewGC.append(idxItem+idx)
            skewCopy.remove(smalestItem)
            idx = idx + 1
        return idxMinSkewGC

    """
        Input   :   1. a genome
                    2. a kmers
                    3. the maximum number of miss match agains the kmers
        Output  :   All the position where kmers appears in genome with at most x missmatches.
    """
    def approxPatternMatching(self,kmers,missMatche):
        idx = 0
        listIdx = []
        globalGenome = self.genome
        while idx <= ( len(globalGenome) - len(kmers)):
            subGenome = globalGenome[idx:(len(kmers)+idx)]
            if subGenome == kmers:
                # We got a match, append this indices to the list.
                # and continuing processing the rest of the genome.
                listIdx.append(idx)
            elif self.isXAtMostMissMatches(subGenome,kmers,missMatche) == 1:
                # checking if this sub-genome has at most x missmatches.
                # YES this subgenome got at most x missmatches, append to the list of indexes.
                listIdx.append(idx)
            idx = idx + 1
        return listIdx

    """
        Input   :   1. a sub-genome
                    2. a kmers
                    3. the max number of miss match between sub-genome and kmers
        Output  :   true if at most x missmatches between sub-genome and kmers
    """
    def isXAtMostMissMatches(self,subGenome,kmers,xMissMatches):
        missMatches = 0
        idx = 0
        for char in subGenome:
            if char == kmers[idx]:
                idx = idx + 1
            else:
                missMatches = missMatches + 1
                if missMatches > xMissMatches:
                    return 0
                idx = idx + 1
        return 1

    def loadCodonTable(self,fileCodonTable):
        """
            Input   :   the file containing the codon table for the genome.
            Output  :   an hash table as key the 3-mers, and as value the amino acids.
        """
        infile  =   open(fileCodonTable,'r')
        listToken   =   []
        for line in infile:
            line = line.replace('\n','')
            listToken = line.split(' ')
            if len(listToken) > 1:
                # We got an amino acids for the codon
                # inserting into the hash table.
                self.codonTable[listToken[0]] =   listToken[1]
            else:
                self.codonTable[listToken[0]] =   '*'

    def translateRNAIntoAcido(self,rna):
        idx = 0
        aminoAcid = ""
        while idx <= ( len(rna) - 3 ):
            codon   =   rna[idx:idx+3]
            acido   =   self.codonTable[codon]
            aminoAcid = aminoAcid + acido
            idx = idx + 3
        return aminoAcid


    def transcriptionDNAToRNA(self):
        """
            Input   :   1. the genome it-self.
            Output  :   The transcription of this genome [ DNA ] into is equivalent RNA.
        """
        rna = ""
        for char in self.genome:
            if char == 'T':
                rna = rna + 'U'
            else:
                rna = rna + char
        return rna