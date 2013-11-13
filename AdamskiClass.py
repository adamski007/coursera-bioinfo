class AdamskiClass:
    """ This is a test for the documentation """
    genome = ""

    def __init__(self, strGenome):
        self.genome = strGenome

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
