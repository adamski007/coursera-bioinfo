class AdamskiClass:
    """ This is a test for the documentation """
    genome = ""

    def __init__(self, strGenome):
        self.genome = strGenome

    def give_complement(self, letter):
        if ( letter == 'A' ):
            return 'T'
        elif ( letter == 'T' ):
            return 'A'
        elif ( letter == 'C' ):
            return 'G'
        elif ( letter == 'G' ):
            return 'C'

    def buildReverseComplement(self):
        complement = ''
        for char in self.genome:
            complement = complement + self.give_complement(char)
        return complement[::-1]

    def findAllOcc(self,kmers):
        listPosition = []
        idx = 0
        while idx <= ( len(self.genome) - len(kmers) ):
            if ( kmers == self.genome[idx:idx+len(kmers)] ) :
                listPosition.append(idx)
            idx = idx + 1
        return listPosition

    def find_x_kmers(self,str_kmers,len_kmers):
        idx = 0;
        count = 0;
        while idx <= ( len(self.genome)-len_kmers):
            if ( self.genome[idx:(idx+len_kmers)] == str_kmers ):
                count = count + 1;
            idx = idx +1;
        return count;

    def findMostFrequentKMers(self,length1):
        # For historic reasons...
        len_kmers = length1
        idx	=	0;
        list_kmers = [];
        most_present_kmers	= 0;
        max_count	= 0;
        while idx <= ( (len(self.genome)) - len_kmers ):
            str_kmers = self.genome[idx:(idx+len_kmers)];
            count = self.find_x_kmers(str_kmers,len_kmers);
            if ( count > max_count ):
                # Re-init the list as we got more a bigger k-mers in this attempt.
                list_kmers = [];
                list_kmers.append(str_kmers);
                max_count = count
            elif ( count == max_count ):
                if not str_kmers in list_kmers:
                    list_kmers.append(str_kmers);
                    # Next line probably not needed. TO DO [ to check actually. ]
                    max_count = count;
            idx = idx + 1;
        return list_kmers