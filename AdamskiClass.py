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
