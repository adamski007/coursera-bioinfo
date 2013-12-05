import sys
import numpy
import os

class AdamskiClass:
    """ This is a test for the documentation """

    def __init__(self, strGenome):
        self.genome = strGenome
        self.codonTable = {}
        self.massTable  =  {}
        self.allMassValue = []
        self.listxmutationkmers = []
        self.overlap_graph = {}
        self.de_bruijn_grapth = {}

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
        all_kmers = self.gene
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

    @staticmethod
    def generateAllKmers(length):
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

    def motifEnumeration(self,listDNA,lengthkmers,xmissmatches):
        """
        Will get an enumeration of the kmers of length lengthkmers
        with at most xmisstaches and present in each dna from listDNA.

        Output : a dictionnary
                    - keys = kmers
                    - values = number of times this kmers was present as a motif.
        """
        allkmersDNA = self.buildingAllKmersFromDna(listDNA,lengthkmers)
        allmotif = {}
        for kmers in allkmersDNA.keys():
            # Iterating over the kmers [ key ] from all dna from listDNA.
            self.generateAllKmersWithXMutation(kmers,xmissmatches)
            allkmersmutated = self.listxmutationkmers[:]
            for kmersmutated in allkmersmutated:
                count = 0
                listIdx = []
                for dna in listDNA:
                    self.genome = dna
                    listIdx = self.approxPatternMatching(kmersmutated,xmissmatches)
                    if len(listIdx) > 0:
                        count+=1
                if count == len(listDNA):
                    # It means that we found the kmers mutated in each dna from listDNA.
                    if allmotif.get(kmersmutated) == None:
                        # the kmers-mutated is not present in the dic, inserting it...
                        allmotif[kmersmutated] = 1
                    else:
                        allmotif[kmersmutated]+=1
        return allmotif


    def buildingAllKmersFromDna(self,listDNA, length):
        """
        As now, we got a DNA from a list, and not anymore from a single line.
        We need a way to get all the possible kmers of a length LENGTH in all dna.
        With a dictionnary, it is faster, because if a kmers already exist, we can find it in O(1)
        the value in each of the key , is the number of kmers we found with his key.
        """
        allkmers = {}
        newkmers = ''
        for dna in listDNA:
            tes = 0
            while idx < len(dna)-length:
                newkmers = dna[idx:idx+length]
                if allkmers.get(newkmers) == None:
                    # kmers not alredy present in the dictionnary, inserting...
                    allkmers[newkmers] = 1
                else:
                    # Adding 1
                    allkmers[newkmers]+=1
                idx+=1
        return allkmers

    @staticmethod
    def rebuildList(listKmers):
        """
        In some cases, we got a list where each element is also a list. This is not quite convenient to paste the result
        back in the field as answer. This function will rebuild as a string each list within this list.
        Finally, we will got a list where each elem is a string of nucleotide.
        This is much better, and easier to paste in the answer field.
        """
        newList = []
        newStr = ''
        for x in listKmers:
            for y in x:
                newStr = newStr + y
            newList.append(newStr)
            newStr = ''
        return newList

    @staticmethod
    def findMotif(pattern,DNA):
        """
            The hamming distance is the number of different nucleotide between pattern and any kmers
            in DNA. pattern is smaller than DNA.
        """
        idx = 0
        motif = 0
        minDiff = len(DNA)
        currDiff = 0
        while idx <= len(DNA)-len(pattern):
            currDiff = AdamskiClass.computeHammingDistance(pattern,DNA[idx:idx+len(pattern)])
            if currDiff < minDiff:
                minDiff = currDiff
                motif = DNA[idx:idx+len(pattern)]
            idx+=1
        return minDiff,motif

    @staticmethod
    def computeHammingDistance(pattern,kmers):
        """
        The hamming distance is the number of differente nucleotide between pattern and a kmers.
        """
        idx = 0
        difference = 0
        while idx < len(pattern):
            if pattern[idx] != kmers[idx]:
                difference+=1
            idx+=1
        return difference

    @staticmethod
    def computeSumHammingDistance(kmers,listDNA):
        """
        For each DNA in listDNA, it will compute the min hamming distance against kmers, and we will
        return the sum of all this.
        """
        distance = 0
        motif =''
        sum = 0
        for dna in listDNA:
            # computing the hamming distance in each of these DNA.
            distance,motif  = AdamskiClass.findMotif(kmers,dna)
            sum = sum + distance
        return sum,motif

    @staticmethod
    def findMedianString(listDNA,lengthKmers):
        """
        The function will find the kmers pattern that minimize distance(pattern,dna) among all
        kmers in each dna of listDNA.
        +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CODE NOT ALWAYS WORKING SUCESSFULLY , TO RE IMPLEMENT AND CHECK.
        +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        """
        allkmers = AdamskiClass.generateAllKmers(lengthKmers)
        # We need to rebuild the list,to get a string in each elem of the list.
        allkmersrebuilded = AdamskiClass.rebuildList(allkmers)
        bestpattern  = ''
        currpattern = ''
        bestscore    = sys.maxsize
        currentscore    =   sys.maxsize
        bestkmers   = ''
        for kmers in allkmersrebuilded:
            currentscore,currpattern = AdamskiClass.computeSumHammingDistance(kmers,listDNA)
            if currentscore < bestscore:
                bestscore   = currentscore
                bestpattern = currpattern
                bestkmers = kmers
        return bestkmers

    @staticmethod
    def readAndBuildListFromFile(nameFile):
        """
        As we now sometime, need to build a list with a DNA in each elem of the list.
        This function will do the trick.
        Building this list manually take too much time, and it is less error prune doing automatically.
        """
        listDNA = []
        infile = open(nameFile,'r')
        for line in infile:
            line = line.replace('\n', '')
            listDNA.append(line)
        return listDNA


    @staticmethod
    def rebuildKmersWithMutation(kmers,newnucleotide,idxnucleotide):
        """
        To re-construct the kmers with his mutated nucleotide,
        we can see 3 different way of doing that.
        # 1. nucleotide mutated is in first place of the kmers.
        # 2. nucleotide mutated is in the middle place of the kmers
        # 3. nucleotide mutated is in the last place of the kmers
        """
        newkmers = ''
        if idxnucleotide == 0:
            newkmers = newkmers + newnucleotide
            newkmers = newkmers + kmers[1:]
        elif idxnucleotide == len(kmers)-1:
            newkmers = newkmers + kmers[:idxnucleotide]
            newkmers = newkmers + newnucleotide
        else:
            # we are in the middle of the kmers
            newkmers = newkmers + kmers[:idxnucleotide]
            newkmers = newkmers + newnucleotide
            newkmers = newkmers + kmers[idxnucleotide+1:]
        return newkmers

    def checkAndInsertNewKmers(self,kmers):
        """
        Function needed with generateAllKmersWithXmutation, with the recursion of the code,
        each time we insert a new kmers, we first check if it is alreaad present, and if not
        we insert the kmers in the list.
        """
        if kmers not in self.listxmutationkmers:
            self.listxmutationkmers.append(kmers)


    def generateAllKmersWithXMutation(self,kmers,xmutation):
        """
        Another methode to try to generate all the kmers that differ from kmers with at most
        x nucleotide mutated.
        """
        newkmers = ''
        for char in 'ATCG':
            idx = 0
            while idx < len(kmers):
                # Mutating a nucleotide in each position of the original kmers.
                newkmers = AdamskiClass.rebuildKmersWithMutation(kmers,char,idx)
                self.checkAndInsertNewKmers(newkmers)
                if xmutation > 1:
                    #print 'doing recursion with kmers : ',newkmers
                    self.generateAllKmersWithXMutation(newkmers,xmutation-1)
                idx+=1

    @staticmethod
    def buildMatrixProfile(matrixfile):
        """
        We will read the file give as input, and build
        his matrix profile. Each line, will represent his vector of proba.
        """
        infile = open(matrixfile)
        listProba = []
        for line in infile:
            line = line.replace('\n','')
            listProba.append(line)
        # We now known the dimensions of the matrix. It will always be 4 column [ A,C,T,G ]
        # and the number of elem in listProba.
        # Initializing the matrix.
        matrixproba = numpy.zeros( (len(listProba),4) )
        # Filling the matrix with elemn in listProba.
        listtoken = []
        idx = 0
        for elem in listProba:
            listtoken = elem.split(' ')
            matrixproba[idx] = listtoken
            idx+=1
        return matrixproba

    @staticmethod
    def getYValueForMatrixProba(nucleotide):
        """
        As in our matrix proba, each column represent a proba for a nucleotide.
        We need to get a way, to get the index according to the nucleotide.
        A -> 0 ; C -> 1 ; G -> 2 ; T -> 3
        The x index is easy, because we parse the file line by line, and hence increase the idx
        x always by one.
        """
        if nucleotide == 'A':
            return 0
        elif nucleotide == 'C':
            return 1
        elif nucleotide == 'G':
            return 2
        elif nucleotide == 'T':
            return 3

    @staticmethod
    def computeProbaKmers(matrixProba, kmers):
        """
        Will compute the proba of this kmers, with the proba of matrixProba
        Example of matrix proba of a kmers of length 3
        A   C   G T
        0   0   0 1     -> proba of 1 for T
        0   1   0 0     -> proba of 1 for C
        0.5 0.5 0 0     -> proba of 0.5 for A , proba of 0.5 for C
        """
        idx = 0
        proba = 1 # as 1 is neutral for a multiplication.
        for nucleotide in kmers:
            x_idx = AdamskiClass.getYValueForMatrixProba(nucleotide)
            proba = proba * matrixProba[idx][x_idx]
            idx+=1
        return proba

    def greedy_Motif_Search(self,listDNA,length_kmers,length_list_dna):
        """
        Slide 39 from chapter 3 of bio-informatics coursera.
        """
        # Building our initial list of kmers for bestMotifs
        list_kmers_best_motifs = []
        for dna in listDNA:
            # Getting the first kmers from each dna of listDNA
            current_kmers = dna[:length_kmers]
            list_kmers_best_motifs.append(current_kmers)
        best_score = AdamskiClass.score_Motifs(list_kmers_best_motifs)
        current_score = sys.maxsize
        idx = 0
        first_dna = listDNA[0]
        motif_one = ''
        list_best_motifs = []
        while idx <= length_kmers:
            current_kmers = first_dna[idx:idx+length_kmers]
            list_motif_for_profile = [current_kmers]
            idx_list_dna = 1
            while idx_list_dna <= length_list_dna-1:
                matrix_count = AdamskiClass.get_Count_Matrix_Motifs(list_motif_for_profile)
                matrix_proba = AdamskiClass.get_Profil_Matrix_Motifs(matrix_count)
                #print idx_list_dna
                #print matrix_proba
                # We need to put the current dna in self.genome, because findMostProbableKmers use that.
                self.genome = listDNA[idx_list_dna]
                """
                print 'The current idx for list dna : ',idx_list_dna
                print 'We will search for the best prob kmers in this dna : ',self.genome
                print 'Matrix profile build with : '
                print list_motif_for_profile
                print '+++++'
                print 'length_kmers : ',length_kmers,'and the matrix proba is : '
                print matrix_proba
                """
                most_probbable_kmers = self.findMostProbableKmers(length_kmers,matrix_proba)
                list_motif_for_profile.append(most_probbable_kmers)
                idx_list_dna+=1
            current_score = AdamskiClass.score_Motifs(list_motif_for_profile)
            #print 'current score : ',current_score, 'and the best score : ',best_score
            #print list_motif_for_profile
            if current_score < best_score:
                list_best_motifs = list_motif_for_profile[:]
                best_score = current_score
            idx+=1
        return list_best_motifs


    def findMostProbableKmers(self,length_kmers,matrixproba):
        """
        The function will the most probably kmers, according to the
        matrix of proba that we already build.
        """
        bestproba = 0
        currentproba = 0
        currentkmers = ''
        # In case, we never find a better kmers, initializing the best kmers
        # as the first we found.
        bestkmers = self.genome[:length_kmers]
        idx = 0
        while idx <= len(self.genome)-length_kmers:
            currentkmers = self.genome[idx:idx+length_kmers]
            currentproba = AdamskiClass.computeProbaKmers(matrixproba, currentkmers)
            if currentproba > bestproba:
                bestkmers = currentkmers
                bestproba = currentproba
            idx+=1
        return bestkmers


    @staticmethod
    def get_Count_Matrix_Motifs(list_kmers):
        """
        Input   : a simple list with all the kmers.
        Output  : a simple matrix computed as describe below.

        It will compute the matrix count of each nucleotide of motifs.
        Motifs        T   C   G   G   G   G   g   T   T   T   t   t
                      c   C   G   G   t   G   A   c   T   T   a   C
                      a   C   G   G   G   G   A   T   T   T   t   C
                      T   t   G   G   G   G   A   c   T   T   t   t
                      a   a   G   G   G   G   A   c   T   T   C   C
                      T   t   G   G   G   G   A   c   T   T   C   C
                      T   C   G   G   G   G   A   T   T   c   a   t
                      T   C   G   G   G   G   A   T   T   c   C   t
                      T   a   G   G   G   G   A   a   c   T   a   C
                      T   C   G   G   G   t   A   T   a   a   C   C

        Count   A:    2   2   0   0   0   0   9   1   1   1   3   0
                C:    1   6   0   0   0   0   0   4   1   2   4   6
                G:    0   0  10  10   9   9   1   0   0   0   0   0
                T:    7   2   0   0   1   1   0   5   8   7   3   4
        """
        # Initializing the matrix, and each kmers in list_kmers should have the same length.
        matrix_count = numpy.zeros( (len(list_kmers[0]),4) )
        counter_nucleotide = {}
        idx_vert = 0
        length_kmers = len(list_kmers[0])
        while idx_vert < length_kmers:
            for kmers in list_kmers:
                current_nucleotide = kmers[idx_vert]
                if counter_nucleotide.get(current_nucleotide) == None:
                    # Inserting the count of this new nucleotide
                    counter_nucleotide[current_nucleotide] = 1
                else:
                    # Nucleotide was already present, adding 1 to his count.
                    counter_nucleotide[current_nucleotide]+=1
                # Inserting now the count of each nucleotide to the matrix count.
                for nucleotide in list(counter_nucleotide.keys()):
                    idx_horiz = AdamskiClass.translate_Nucleotide_To_Horiz_Idx(nucleotide)
                    #matrix_count[idx_horiz][idx_vert] = counter_nucleotide[nucleotide]
                    matrix_count[idx_vert][idx_horiz] = counter_nucleotide[nucleotide]
                # We should clear the dict for the next cycle.
            counter_nucleotide.clear()
            idx_vert+=1
        return matrix_count

    @staticmethod
    def get_Profil_Matrix_Motifs(matrix_count):
        """
        It will compute the profile matrix with the frequency count of the matrix received
        as input.
        Input   :   a matrix with frequency count for each nucleotide
        Output  :   a matrix where each elem is a proba of this count, relative to each column.
        Count Matrix    A:   2   2   0   0   0   0   9   1   1   1   3   0
                        C:   1   6   0   0   0   0   0   4   1   2   4   6
                        G:   0   0  10  10   9   9   1   0   0   0   0   0
                        T:   7   2   0   0   1   1   0   5   8   7   3   4

        Profile Matrix  A:  .2  .2   0   0   0   0  .9  .1  .1  .1  .3   0
                        C:  .1  .6   0   0   0   0   0  .4  .1  .2  .4  .6
                        G:   0   0   1   1  .9  .9  .1   0   0   0   0   0
                        T:  .7  .2   0   0  .1  .1   0  .5  .8  .7  .3  .4
        """
        # We will count one time the total nucleotide per column,
        # as each time it will be the same.
        total_count = 0
        idx = 0
        # We known in advance that the matrix is max 4 column, because of ACGT
        while idx < 4:
            total_count+=matrix_count[0][idx]
            idx+=1
        # Normalizing the matrix to get proba instead of frequency count.
        return matrix_count/total_count

    @staticmethod
    def translate_Nucleotide_To_Horiz_Idx(nucleotide):
        """
        In the matrix count, or profile matrix, the horizontal idx are :
        A -> 0
        C -> 1
        G -> 2
        T -> 3
        and the vertical idx are increased as we parse the kmers.
        """
        if nucleotide == 'A':
            return 0
        elif nucleotide == 'C':
            return 1
        elif nucleotide == 'G':
            return 2
        elif nucleotide == 'T':
            return 3

    @staticmethod
    def score_Motifs(list_kmers):
        """
        If kmers are put in a line, a score of a list of kmers is defined as :
        the sum of the number of un-popular nucleotide in each column.
        Motifs   [    T   C   G   G   G   G   g   T   T   T   t   t
                      c   C   G   G   t   G   A   c   T   T   a   C
                      a   C   G   G   G   G   A   T   T   T   t   C
                      T   t   G   G   G   G   A   c   T   T   t   t
                      a   a   G   G   G   G   A   c   T   T   C   C
                      T   t   G   G   G   G   A   c   T   T   C   C
                      T   C   G   G   G   G   A   T   T   c   a   t
                      T   C   G   G   G   G   A   T   T   c   C   t
                      T   a   G   G   G   G   A   a   c   T   a   C
                      T   C   G   G   G   t   A   T   a   a   C   C     ]

        Score         3 + 4 + 0 + 0 + 1 + 1 + 1 + 5 + 2 + 3 + 6 + 4 = 30
        The best score, would be the list of kmers that minimize this sum.
        PS : There is the same number of nucleotide in each kmers of the list.
        """
        length_kmers = len(list_kmers[0])
        idx = 0
        counter_nucleotide = {}
        list_unpopular_count = []
        while idx < length_kmers:
            # we will scan each kmers, and compute the popular nucleotide in each column.
            for kmers in list_kmers:
                current_nucleotide = kmers[idx]
                if counter_nucleotide.get(current_nucleotide) == None:
                    # Inserting the count of this new nucleotide
                    counter_nucleotide[current_nucleotide] = 1
                else:
                    # Nucleotide was already present, adding a 1 to his count.
                    counter_nucleotide[current_nucleotide]+=1
            # Adding the count of un-popular nucleotide in this column in a list
            count_unpopular_nucleotide = AdamskiClass.counting_UnPopular_Nucleotide(counter_nucleotide)
            list_unpopular_count.append(count_unpopular_nucleotide)
            # Clear the dict for the next column
            counter_nucleotide.clear()
            idx+=1
        # Adding all the un-popular count in the list to get the score of the motifs.
        sum = 0
        for unpopular_count in list_unpopular_count:
            sum+=unpopular_count
        return sum


    @staticmethod
    def counting_UnPopular_Nucleotide(counter_nucleotide):
        """
        Method used for scoring a list of motif.
        In the dictionnary [ counter_nucleotide ], there is a count for each
        nucleotide, and the method will return the number of un-popular nucleotide.
        """
        total_count_nucleotide = 0
        popular_nucleotide_count = 0
        popular_nucleotide = ''
        # 1. Find the most popular nucleotide.
        for key in list(counter_nucleotide.keys()):
            total_count_nucleotide+=counter_nucleotide[key]
            if counter_nucleotide[key] > popular_nucleotide_count:
                # We got the most popular nucleotide.
                popular_nucleotide_count = counter_nucleotide[key]
                popular_nucleotide = key
        # 2. Getting the number of un-popular nucleotide to get it's score.
        return total_count_nucleotide - popular_nucleotide_count



    def findStringComposition(self,lengthkmers):
        """
        This function will output [ in a list ] all the kmers [ sub-string ] of DNA,
        even the one repeated, in a lexicographic order.
        """
        kmerssorted = []
        idx = 0
        newkmers = ''
        while idx <= len(self.genome)-lengthkmers:
            newkmers = self.genome[idx:idx+lengthkmers]
            kmerssorted.append(newkmers)
            idx+=1
        kmerssorted.sort()
        return kmerssorted

    @staticmethod
    def get_prefixe_kmers(kmers):
        """
        The prefixe of a kmers is all it's character except the last one.
        """
        return kmers[:len(kmers)-1]

    @staticmethod
    def get_suffix_kmers(kmers):
        """
        The suffixe of a kmers is all it's character except the first one.
        """
        return kmers[1:]


    def build_Overlap_Graph(self,list_kmers):
        """
        With all the kmers present in list_kmers , we will construct a dictionnary where :
            key     : is a kmers
            value   : another kmers who prefixe is the suffixe of the kmers as the key
        By instance :
        ATGCG       AGGCA
        GCATG       GGCAT
        CATGC
        Dictionnary key         dictionnary value
        AGGCA           ->      GGCAT
        CATGC           ->      ATGCG
        GCATG           ->      CATGC
        GGCAT           ->      GCATG
        """
        # For each kmers_suf in list_kmers , we need to find another kmers_pre
        # which got Suffixe( kmers_suf ) = Prefixe( kmers_pre )
        # and add that to the dictionnary, if not already present.
        for kmers_suf in list_kmers:
            for kmers_pre in list_kmers:
                if self.get_suffix_kmers(kmers_suf) == self.get_prefixe_kmers(kmers_pre):
                    # Checking and inserting.
                    if self.overlap_graph.get(kmers_suf) == None:
                        # Value does not exist, inserting...
                        self.overlap_graph[kmers_suf] = [kmers_pre]
                    else:
                        # Value does already exist, so getting the present value
                        new_list = self.overlap_graph[kmers_suf]
                        new_list.append(kmers_pre)
                        self.overlap_graph[kmers_suf] = new_list

    def print_edge_overlap_graph(self):
        """
        This is the function used to print the overlap grapth
        the way bio-informatics courses want us to do.
        KMERS_1 -> KMERS_2
        """
        # We first need all the key of the dic in a list, this way we will be able
        # to sort that list, and print it in a lexicograph way. Like string composition.
        all_keys = list(self.overlap_graph.keys())
        all_keys.sort()
        for key in all_keys:
            all_kmers_values_list = self.overlap_graph[key]
            all_kmers_values_str  = ''
            for kmers in all_kmers_values_list:
                #all_kmers_values_str = kmers + ' ' + all_kmers_values_str
                # the white-space at the end , give error in stepic answer ?
                all_kmers_values_str = kmers
            # We now got all the value of a particular keys in a simple string
            print key,'->',all_kmers_values_str

    def build_DeBruin_Graph(self,stringTXT,length_kmers):
        """
        Slide 53 of Chapter 4.
        This is the same kind of implementation of the overlapping graph, except
        that here, we got as input a string and a length of kmers.
        So :
        1. we should decompose the string in slices kmers of k length.
        2. Build the overlap-graph, as it takes already into account the fact
            that more than 1 kmers as prefixe of a suffixe.
        3. From the graph at point two, re-build the same but remove always the
            last char of each node.
            A node is a pair : key / value
                                where value is a list of kmers.
            So, remove the last char from key, and from each elem in the list value.
        Sample Input:
             4
             AAGATTCTCTAC

        Sample Output:
             AAG -> AGA
             AGA -> GAT
             ATT -> TTC
             CTA -> TAC
             CTC -> TCT
             GAT -> ATT
             TCT -> CTA,CTC
             TTC -> TCT
        """
        self.genome = stringTXT
        # Point 1 execution.
        dna_composition_into_kmers = self.findStringComposition(length_kmers)
        # Building the overlap-grapth.
        self.build_Overlap_Graph(dna_composition_into_kmers)
        # The graph is now stored in [ self.overlap_graph ]
        # Point 3, iterating over overlap-graph, and removing each time last nucleotide.
        for key in list(self.overlap_graph.keys()):
            key_shortened = key[:len(key)-1]
            list_value_kmers_from_key_shortened = []
            list_kmers = list(self.overlap_graph[key])
            for kmers in list_kmers:
                kmers_shortened = kmers[:len(kmers)-1]
                list_value_kmers_from_key_shortened.append(kmers_shortened)
            self.de_bruijn_grapth[key_shortened] = list_value_kmers_from_key_shortened

    def print_edge_debruijn_graph(self):
        """
        Again the course want us to print the debruijn graph in a
        lexicograph way.
        So, we need to sort the keys, and the print their related content.
        """
        keys_sorted = list(self.de_bruijn_grapth.keys())
        keys_sorted.sort()
        output_string = ''
        for key in keys_sorted:
            output_string = key + ' -> '
            all_kmers_associated_with_key = list(self.de_bruijn_grapth[key])
            idx = 0
            # We will parse all_mers_associated_with_key with an index,
            # this way we known when the list is finished, because we don t
            # have to add an [ , ] at the end of the list, only between
            # 2 kmers.
            while idx < len(all_kmers_associated_with_key):
                kmers = all_kmers_associated_with_key[idx]
                if idx == len(all_kmers_associated_with_key)-1:
                    # This is the last element of the list
                    # Don't add a [ , ] to the print, just the last elem
                    output_string = output_string + kmers
                else:
                    # We need to add [ , ] to the oupput, as there is
                    # still some elem to be printed.
                    output_string = output_string + kmers + ','
                idx+=1
            print output_string

    def build_List_Decomposition_DNA_Sequential(self,lengthkmers):
        """
        It will decompose DNA into kmers
        # TO FINISHE YET !!!!!!!!!!!!!!!!!!!!!
        """
        list_kmers = []
        idx = 0
        newkmers = ''
        while idx <= len(self.genome)-lengthkmers:
            newkmers = self.genome[idx:idx+lengthkmers]
            list_kmers.append(newkmers)
            idx+=1
        return list_kmers

    def build_DeBruijn_Grapth(self,length_kmers,dna):
        """
        Will build the algo mentionned on slide 53 chapter 4
        """

















