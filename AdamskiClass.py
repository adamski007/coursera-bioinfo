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
