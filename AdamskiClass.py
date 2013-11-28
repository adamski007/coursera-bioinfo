class AdamskiClass:
    """ This is a test for the documentation """

    def __init__(self, strGenome):
        self.genome = strGenome
        self.codonTable = {}
        self.massTable  =  {}
        self.allMassValue = []

    def buildAllMassValue(self):
        """
        Building a list with all the mass possible.
        Easier to check if a mass exist.
        """
        for mass in self.massTable:
            # As all mass are an integer, we convert directly to integer, otherwise python think that it is a string.
            self.allMassValue.append(int(self.massTable[mass]))


    def give_complement(self, letter):
        """
        Input : 1. a letter in the genome.
        Output: The complement of this letter in the genome.
        """
        if ( letter == 'A' ):
            return 'T'
        elif ( letter == 'T' ):
            return 'A'
        elif ( letter == 'C' ):
            return 'G'
        elif ( letter == 'G' ):
            return 'C'


    def buildReverseComplement(self):
        """
            Input : 1. a genome
            Ouput : the reverse complement of this genome.
        """
        complement = ''
        for char in self.genome:
            complement = complement + self.give_complement(char)
        return complement[::-1]


    def findAllOcc(self, kmers):
        """
            Input : 1. a kmers
            Output: a list with the position where the kmers a located in the genome.
        """
        listPosition = []
        idx = 0
        while idx <= ( len(self.genome) - len(kmers) ):
            if ( kmers == self.genome[idx:idx + len(kmers)] ):
                listPosition.append(idx)
            idx = idx + 1
        return listPosition

    def find_x_kmers(self, str_kmers, len_kmers):
        """ Input :     1. a genome
                        2. a kmers
            Output:     the numbers of kmers we found in the genome.
        """
        idx = 0;
        count = 0;
        while idx <= ( len(self.genome) - len_kmers):
            if ( self.genome[idx:(idx + len_kmers)] == str_kmers ):
                count = count + 1;
            idx = idx + 1;
        return count;

    def findMostFrequentKMers(self, length1):
        """
            Input   :   1. a genome
                        2. length of a kmers
            Output  :   a list containing all the kmers with a specific size in the genome.
        """
        # For historic reasons...
        len_kmers = length1
        idx = 0
        list_kmers = []
        most_present_kmers = 0
        max_count = 0
        while idx <= ( (len(self.genome)) - len_kmers ):
            str_kmers = self.genome[idx:(idx + len_kmers)]
            count = self.find_x_kmers(str_kmers, len_kmers)
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


    def findClump(self, size_window, size_kmers, number):
        """
            Input   :   1. a genome
                        2. a size of window where to find the number of kmers
                        3. the length of the kmers
                        4. the number of kmers we need to find the specified window.
            Output  :   The list of kmers found in the window.
        """
        idx = 0
        list_kmers = []
        globalGenome = self.genome
        while idx <= ( len(globalGenome) - size_window ):
            window = self.genome[idx:(idx + size_window)]
            self.genome = window
            # Doing the computation temporarily
            list_kmers_in_window = self.findMostFrequentKMers(size_kmers)
            # Checking if these kmers are present enough in the window.
            for kmers in list_kmers_in_window:
                if self.find_x_kmers(kmers, len(kmers)) >= number:
                    # It does means that we found at least the requested number of kmers in
                    # the specified window. We got one candidate for the the kmers.
                    if not kmers in list_kmers:
                        list_kmers.append(kmers)
            idx = idx + 1
            self.genome = globalGenome
        return list_kmers


    def computeSkewGC(self):
        """
            Input   : 1. a genome
            Output  : a list containing at each idx of the genome, the current skew
                The skew is defined as the defirence between G and C at each indices of the genome.
        """
        skewGC = [0]
        curSkew = 0
        for char in self.genome:
            if char == 'C':
                curSkew = skewGC[-1] - 1
                skewGC.append(curSkew)
            elif char == 'G':
                skewGC.append(skewGC[-1] + 1)
            else:
                skewGC.append(skewGC[-1])
        return skewGC

    def getIdxMinSkewGC(self, listSkew):
        """
            Input   : the skew of a genome
            Output  : a list containing the index where the skew is the minimum.
        """
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
            idxMinSkewGC.append(idxItem + idx)
            skewCopy.remove(smalestItem)
            idx = idx + 1
        return idxMinSkewGC

    def approxPatternMatching(self, kmers, missMatche):
        """
            Input   :   1. a genome
                        2. a kmers
                        3. the maximum number of miss match agains the kmers
            Output  :   All the position where kmers appears in genome with at most x missmatches.
        """
        idx = 0
        listIdx = []
        globalGenome = self.genome
        while idx <= ( len(globalGenome) - len(kmers)):
            subGenome = globalGenome[idx:(len(kmers) + idx)]
            if subGenome == kmers:
                # We got a match, append this indices to the list.
                # and continuing processing the rest of the genome.
                listIdx.append(idx)
            elif self.isXAtMostMissMatches(subGenome, kmers, missMatche) == 1:
                # checking if this sub-genome has at most x missmatches.
                # YES this subgenome got at most x missmatches, append to the list of indexes.
                listIdx.append(idx)
            idx = idx + 1
        return listIdx


    def mostFrequentKmersXMissMatche(self,len_kmers,xmissMatches):
        """
            Will find the most present kmers present in the DNA, with at most x miss-matches.
            Input   :   1. The genome it-self
                        2. the choosen len of a kmers
                        3. the max number of miss-matches agains the kmers.
            Output  :   The most frequent kmers with at most x miss-matches.
        """
        idx = 0
        countKMersPresent = 0
        listKmersMostPresent = []
        while idx <= len(self.genome)-len_kmers:
            currentKMers    = self.genome[idx:idx+len_kmers]
            listIdx = self.approxPatternMatching(currentKMers,xmissMatches)
            print currentKMers,len(listIdx),listIdx
            if len(listIdx) > countKMersPresent:
                # We got a new record for the most present kmers.
                listKmersMostPresent = []
                #kmers = self.genome[listIdx[0]:listIdx[0]+len_kmers]
                listKmersMostPresent.append(currentKMers)
                countKMersPresent   = len(listIdx)
            elif len(listIdx) == countKMersPresent and len(listIdx) > 0:
                # We got a [new] kmers
                kmers = self.genome[listIdx[0]:listIdx[0]+len_kmers]
                if currentKMers not in listKmersMostPresent:
                    listKmersMostPresent.append(currentKMers)
            idx = idx + 1
        return listKmersMostPresent

    def mostFrequentKmersXMissMatches2(self,len_kmers,xMissMatches):
        """
            Will find the most present kmers present in the DNA, with at most x miss-matches.
            Input   :   1. The genome it-self
                        2. the choosen len of a kmers
                        3. the max number of miss-matches agains the kmers.
            Output  :   The most frequent kmers with at most x miss-matches.
        """
        # This version will take into account the kmers that even not present in the DNA it-self.
        print 'Generating all kmers of length : ',len_kmers
        all_kmers = self.generateAllKmers(len_kmers)
        print all_kmers[0]
        print 'Generation of all the kmers done.'
        countKMersPresent = 0
        listKmersMostPresent = []
        print len(all_kmers)
        for kmers in all_kmers:
            #print kmers
            listIdx = self.approxPatternMatching(kmers,xMissMatches)
            #print listIdx
            if len(listIdx) > countKMersPresent:
                # We got a new record for the most present kmers.
                listKmersMostPresent = []
                #kmers = self.genome[listIdx[0]:listIdx[0]+len_kmers]
                listKmersMostPresent.append(kmers)
                countKMersPresent   = len(listIdx)
            elif len(listIdx) == countKMersPresent and len(listIdx) > 0:
                # We got a [new] kmers
                if kmers not in listKmersMostPresent:
                    listKmersMostPresent.append(kmers)
        return listKmersMostPresent


    def isXAtMostMissMatches(self, subGenome, kmers, xMissMatches):
        """
            Input   :   1. a sub-genome
                        2. a kmers
                        3. the max number of miss match between sub-genome and kmers
            Output  :   true if at most x missmatches between sub-genome and kmers
        """
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


    def loadCodonTable(self, fileCodonTable):
        infile = open(fileCodonTable, 'r')
        listToken = []
        for line in infile:
            line = line.replace('\n', '')
            listToken = line.split(' ')
            if len(listToken) > 1:
                # We got an amino acids for the codon
                # inserting into the hash table.
                self.codonTable[listToken[0]] = listToken[1]
            else:
                self.codonTable[listToken[0]] = '*'


    def translateRNAIntoAcido(self, rna):
        idx = 0
        aminoAcid = ""
        while idx <= ( len(rna) - 3 ):
            codon = rna[idx:idx + 3]
            acido = self.codonTable[codon]
            aminoAcid = aminoAcid + acido
            idx = idx + 3
        return aminoAcid


    def transcribeDNAToRNA(self):
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


    def checkSubStr(self, subStr, acido):
        globalGenome = self.genome
        self.genome = subStr
        rna = self.transcribeDNAToRNA()
        rnaComp = self.buildReverseComplement()
        self.genome =   rnaComp
        rnaRevComp = self.transcribeDNAToRNA()
        self.genome = globalGenome
        if self.translateRNAIntoAcido(rna) == acido:
            # we got a match
            return 1
        elif self.translateRNAIntoAcido(rnaRevComp) == acido:
            # checking with the complement
            return 1
        else:
            return 0


    def findSubStringAcido(self, acido):
        """
                Input   :   1. the genome
                            2. the acido which was encoded from a sub-string of the genome.
                Output  :   a list of idx of all sub-string of the genome which encode a particular acido amids.
        """
        idx = 0
        listMatches = []
        while idx <= len(self.genome) - len(acido) * 3:
            subStr = self.genome[idx:idx + len(acido) * 3]
            if 1 == self.checkSubStr(subStr, acido):
                listMatches.append(subStr)
            idx = idx + 1
        return listMatches

    def generateListSubPep(self,amino):
        """
            Out comment here...
            the pseudo algo should be something like that...
            1. from 0 -> len , find the spectruym
               len of 2 , for the last one idx[len-1,len] + len[0:1]
               from 3 -> x :

            2. and for each sub-peptide of len x , give his spectrum, to be done in another function....
        """
        # This list will contains all the subpeptide of amino.
        listSubPep = []
        # Adding the acido amino it-self, and all it's char.
        # Then, we will need to add all sub-pep possible.
        listSubPep.append(amino)
        for char in amino:
            listSubPep.append(char)
        # max len of any subpeptide.
        maxLenSubPep    =   len(amino)-1
        # we must begin with a lenth of subpeptide of min 1.
        curLenSubPep    =   2
        idx     =   0
        while curLenSubPep <= maxLenSubPep:
            # Doing all the length possible for sub-pep.
            while idx < len(amino):
                if idx > len(amino)-curLenSubPep:
                    # do special processing
                    remain = len(amino) - idx
                    subPepRear = amino[idx:idx+remain]
                    addFront = curLenSubPep - remain
                    subPepFront = amino[0:addFront]
                    subPepAll = subPepRear + subPepFront
                    listSubPep.append(subPepAll)
                else:
                    listSubPep.append(amino[idx:idx+curLenSubPep])
                idx =   idx + 1
            curLenSubPep = curLenSubPep + 1
            idx = 0
        return listSubPep

    def computeMassSpectrum(self,listSubPep):
        listMassSpectrum = []
        currentMass = 0
        # We always need to add a mass of zero in the spectrum.
        listMassSpectrum.append(currentMass)
        for subPep in listSubPep:
            # compute the mass spectrum of each elem of the list.
            for char in subPep:
                massElem    =   self.massTable[char]
                currentMass =   currentMass + int(massElem)
            listMassSpectrum.append(currentMass)
            currentMass = 0
        listMassSpectrum.sort()
        return listMassSpectrum

    def loadTableSpectrum(self,fileSpectrum):
        """
            It will load the file [ fileSpectrum ] as the basic spectrum needed
            to generate the spectrum of any amino acid.
            Input   :   the file containing the table of each basic mass spectrum
            Output  :   the hashtable containng all mass spectrum
        """
        infile = open(fileSpectrum, 'r')
        listToken = []
        for line in infile:
            line = line.replace('\n', '')
            listToken = line.split(' ')
            if len(listToken) > 1:
                # We got an amino acids for the codon
                # inserting into the hash table.
                self.massTable[listToken[0]] = listToken[1]
            else:
                self.massTable[listToken[0]] = '*'

    def computeSpectralConvolution(self,spectrum):
        """
        We build a matrix, where each point i , j is a difference between two masses of a spectrum.
        Intput  :   a mass spectrum
        Output  :   a list of element between 57 and 200 DA.
        """
        idxHoriz=   0
        # On the horizontale idx, the last value [ the complete mass of the peptide ] is not present.
        idxVert =   1 + idxHoriz
        # We don't need the elem of the diagonal of the matrix, that's why idxVert begins with 1.
        # As for the computation in the matrix on the right side, we always need less point, we increase the
        # idxVert by one each time.
        listElem    =   []
        while idxHoriz <= len(spectrum)-2:
            while idxVert <= len(spectrum)-1:
                difference  =   int(spectrum[idxVert])-int(spectrum[idxHoriz])
                if difference > 0:
                    # A spectrum of zero is not relevant, that is why we add only the element where is above 0
                    listElem.append(difference)
                # computing next elem verticaly.
                idxVert =   idxVert + 1
            # computing next elem horizontaly
            idxHoriz    =   idxHoriz + 1
            # re-setting the vertical indexe for the next iteration, and increasing the value by idxHoriz
            # because next time we compute [ next column ] , the computation occured one element lower in the matrix,
            # do a picture and you will understand...
            # This is because we compute only the elements under the diagonale of the matrix.
            idxVert =   1 + idxHoriz
        return listElem

    def cyclopeptideSequencing(self,spectrum):
        firstList   =   self.findOneMers(spectrum)
        firstList.append(0)
        while len(firstList) > 0:
            print 'toto'
        return 0

    def findOneMers(self,spectrum):
        """
        Find all the 1-mers compatible with the spectrum
        """
        listOneMers = []
        for mass in spectrum:
            if mass in self.allMassValue:
                # this is a starting point...
                listOneMers.append([mass])
        return listOneMers

    def expandList(self,listMass):
        newListMass = []
        for elem in listMass:
            for mass in self.allMassValue:
                # extending all peptide with each acido.
                tmpList =   elem.append(mass)
                newListMass.append(tmpList)
        return newListMass

    def generateAllKmers(self,length):
        """
        Will generate a list containing all the kmers of length LENGTH.
        Needed for one exercice in coursera.
        Input   :   the length wanted of a kmers
        Output  :   a list with all kmers possible of this length.
        """
        idx = 1
        listKmers = []
        listKmers2 = []
        # Initializing the list with one letter.
        for char in 'ATCG':
            listKmers.append([char])
        idx = idx + 1
        while idx <= length:
            for list in listKmers:
                for char in 'ATCG':
                    newList = list + [char]
                    listKmers2.append(newList)
            listKmers = listKmers2
            listKmers2 = []
            idx = idx + 1
        return listKmers































