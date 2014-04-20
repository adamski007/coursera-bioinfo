import sys
import datetime
import time
import itertools
import pprint
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
        self.count_in_out_edges = {}
        self.edges_weigth_dag = {}
        self.predecessor_nodes = {}
        self.weight_nodes = {}
        self.graph_dag = {}
        self.node_used = []
        self.idx_matrix_amino_acid = {}
        # In the matrix, at [0,0] , this is always a minimal score.
        self.max_value_matrix_score = (0, 0)
        # Initializing the matrix, but will be updated with the current request made on demand.
        self.matrix_multiple_alignment = numpy.zeros((3, 3))
        self.matrix_backtracking_multiple = numpy.zeros((3, 3))
        # Variables needed for the tries construction
        self.index_tries_construction = 1
        self.tries_construction = {}
        self.deep_level_in_tries = 0
        self.deepest_level_longest_repeat = 0
        self.longest_repeat = ''
        # END : Variables needed for the tries construction
        # Variables needed for the suffix tree structure.
        self.num_letter_down_tree = 0
        # END : Variables needed for the suffix tree structure
        # Variables needed for the BW transform.
        self.idx_character_last_column = {}
        # END : Variables needed for the BW transform.
        # Variables needed for better BW_matching .
        self.bw_count = numpy.zeros( (1,1) )
        self.count_char = [0,0,0,0,0]

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


    def findClump(self, size_kmers, size_window, number):
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

    @staticmethod
    def nucleotide_count(dna_string):
        """
            The function will output the count of each nucleotide [ A , C , G , T ] present in a DNA string.
        """
        # First element in the following list will count the number of nucleotide for A, second for C , ...
        # in the DNA string.
        nucleotide_count = [0,0,0,0]
        # We are only working with capital character...
        dna_string_upper = dna_string.upper()
        for nucleotide in dna_string_upper:
            if nucleotide == 'A':
                nucleotide_count[0] = nucleotide_count[0] + 1
            elif nucleotide == 'C':
                nucleotide_count[1] = nucleotide_count[1] + 1
            elif nucleotide == 'G':
                nucleotide_count[2] = nucleotide_count[2] + 1
            elif nucleotide == 'T':
                nucleotide_count[3] = nucleotide_count[3] + 1
            else:
                print("We should not see this lines -> ERROR")
        return nucleotide_count


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
            print(currentKMers,len(listIdx),listIdx)
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
        print(all_kmers[0])
        print('Generation of all the kmers done.')
        countKMersPresent = 0
        listKmersMostPresent = []
        print(len(all_kmers))
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
        return self.genome.replace('T','U')


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

    def computeMassSpectrum(self, listsubpep):
        """
            From the list as input [ listsubpep ], we will build for each of them [ elem of the list ]
            their respectif mass according to the masstable monoisotopic.
            The function will return as output a list, where each element his the mass of each
            corresponding element in the list listsubpep.
            Example :
            Input   : [ SKADYEK ]
            Output  : [ 821.392 ]
            The monoisotopic masstable can be found on wikipedia :
            - http://en.wikipedia.org/wiki/Proteinogenic_amino_acid#Mass_spectrometry
        """
        listMassSpectrum = []
        currentMass = 0.0
        # We always need to add a mass of zero in the spectrum.
        listMassSpectrum.append(currentMass)
        for subPep in listsubpep:
            # compute the mass spectrum of each elem of the list.
            for char in subPep:
                massElem    =   self.massTable[char]
                currentMass =   currentMass + float(massElem)
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
            listToken = line.split()
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
            print('toto')
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
    def compute_fibonacci_rabbit(month_needed, multiple_rabbit, list_relation_recurrence):
        """
            Fn  = Fn-1 + Fn-2
            The number of rabbit in generation n is the number of rabbit at generation n-1
            plus a multiple of the rabbit present in the generation -2.
            Example :
            +++++++++
            F1	F2	F3	F4	F5	F6	F7	F8
            1	1	2	2	4	6	10	16      where n = 8 and k = 1
            1	1	4	7	19	40	97	217     where n = 8 and k = 3
            Example call function : compute_fibonacci_rabbit(5,3,[1,1]) will produce : 19
        """
        if len(list_relation_recurrence) == month_needed:
            return list_relation_recurrence[-1]
        else:
            # Computing the next month number of rabbit needed.
            new_rabbit_count = list_relation_recurrence[-1] + multiple_rabbit * list_relation_recurrence[-2]
            list_relation_recurrence.append(new_rabbit_count)
            return AdamskiClass.compute_fibonacci_rabbit(month_needed, multiple_rabbit, list_relation_recurrence)


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
    def find_motif(self, dna='', kmers=''):
        """
            The function will return a list where each element is the index where we can find [ kmers ]
            in the dna it-self.
            Example :
            Input   : ACTTGCT , CT
            Output  : [ 2 , 6 ]
        """
        list_idx_kmers = []
        idx = dna.find(kmers) + 1  # +1 because starting number in a string is 0 in python.
        if idx != 0:
            list_idx_kmers.append(idx)
            list_sub_idx = AdamskiClass.find_motif(dna[idx+1:], kmers)
            if len(list_sub_idx) == 0:
                return list_idx_kmers
            else:
                list_idx_kmers.extend(list_sub_idx)
                return list_idx_kmers
        else:
            return []


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
    def readAndBuildListPermutation(nameFile):
        """
        Function to read the content of a file which contains a list of signed number, and
        build a list with these numbers.
        """
        list_number = []
        list_number_int = []
        infile = open(nameFile,'r')
        for line in infile:
            line = line.replace('\n', '')
            list_number = line.split(' ')
            for x in list_number:
                list_number_int.append(int(x))
        return list_number_int

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
    def get_consensus_profile( dic_fasta):
        """
        All the dna will be stocked in a FASTA file. The file has already been parsed, and put
        in a dictionnary. What we need to do is to re-execute the method above :
        [ get_Count_Matrix_Motifs ] and transpose the matrix, to get it the way the instructor want.
        And then checking which index got the max number, and according to that printing the proper
        nucleotide.
        0 -> A ; 1 -> C ; 2 -> G ; 3 -> T.
        One sample dataset [ FASTA FILE ] :
                >Rosalind_1
                ATCCAGCT
                >Rosalind_2
                GGGCAACT
                >Rosalind_3
                ATGGATCT
                >Rosalind_4
                AAGCAACC
                >Rosalind_5
                TTGGAACT
                >Rosalind_6
                ATGCCATT
                >Rosalind_7
                ATGGCACT
        will output :
                ATGCAACT
                A: 5 1 0 0 5 5 0 0
                C: 0 0 1 4 2 0 6 1
                G: 1 1 6 3 0 1 0 0
                T: 1 5 0 0 0 1 1 6
        """
        matrix_count = AdamskiClass.get_Count_Matrix_Motifs( dic_fasta.values())
        matrix_transposed = matrix_count.transpose()
        idx_dna = matrix_transposed.argmax(0)
        consensus_dna = ''
        for idx in idx_dna:
            if idx == 0:
                consensus_dna = consensus_dna + 'A'
            elif idx == 1:
                consensus_dna = consensus_dna + 'C'
            elif idx == 2:
                consensus_dna = consensus_dna + 'G'
            elif idx == 3:
                consensus_dna = consensus_dna + 'T'
        print(consensus_dna)
        print('A: ', matrix_transposed[0])
        print('C: ', matrix_transposed[1])
        print('G: ', matrix_transposed[2])
        print('T: ', matrix_transposed[3])


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

    def build_overlap_grapth_from_file(self,nameFile):
        """
        Getting the value of the overlap from a file as :
        0001 -> 0011
        0011 -> 0110
        ...
        """
        self.overlap_graph.clear()
        infile = open(nameFile,'r')
        for line in infile:
            line = line.replace('\n', '')
            list_token = line.split(' -> ')
            key = list_token[0]
            all_values = list_token[1].split(',')
            for value in all_values:
                # Checking and inserting.
                if self.overlap_graph.get(key) == None:
                    # Value does not exist, inserting...
                    self.overlap_graph[key] = [value]
                else:
                    # Value does already exist, so getting the present value
                    new_list = self.overlap_graph[key]
                    new_list.append(value)
                    self.overlap_graph[key] = new_lis


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
            print(key,'->',all_kmers_values_str)

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
            print(output_string)

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
        self.genome = dna
        list_kmers_sequential = []
        list_kmers_sequential = self.build_List_Decomposition_DNA_Sequential(length_kmers)
        idx = 0
        # until len(list_kmers)-2 because the last element should be the value,
        # and the before last the key
        while idx <= len(list_kmers_sequential)-2:
            key = list_kmers_sequential[idx]
            value = list_kmers_sequential[idx+1]
            if self.de_bruijn_grapth.get(key) == None:
                # key still does not exist, creating and inserting the value.
                self.de_bruijn_grapth[key] = [value]
            else:
                # key already exist, adding a value to the list
                cur_list = self.de_bruijn_grapth[key]
                cur_list.append(value)
                self.de_bruijn_grapth[key] = cur_list
            idx+=1

    def test_Insert_kmers_Into_DeBruijn_Graph(self,key,value):
        if self.de_bruijn_grapth.get(key) == None:
            # key still does not exist, creating and inserting the value.
            self.de_bruijn_grapth[key] = [value]
        else:
            # key already exist, adding a value to the list
            cur_list = self.de_bruijn_grapth[key]
            cur_list.append(value)
            self.de_bruijn_grapth[key] = cur_list

    def build_deBruijn_graph_paired_reads(self):
        """
        TO DO - TO IMPLEMENT STILL
        READ THE FILE FIRST WITH THE RELEVANT DATA.
        """

    def build_DeBruijn_Graph_from_listKmers(self,list_kmers,namefile=''):
        """
        The goal is also here to build this De Bruijn graph, but from
        a list of kmers.
        The task is here to build an overlap graph from these kmers,
        and from each key, value build the edges.
        By instance, for the overlap node :
        - ATG -> TGC
        We got finally :
        - AT -> TG
        - TG -> GC
        """
        # Building our overlap graph, and from that building debruijn graph.
        start_time = datetime.datetime.now()
        if namefile != '':
            # We got the overlap graph in a file, don t needed to be build it.
            self.build_overlap_grapth_from_file(namefile)
        else:
            self.build_Overlap_Graph(list_kmers)
        end_time   = datetime.datetime.now()
        print('Overlap graph has been build in                      : ',end_time-start_time)
        list_kmers_processed = []
        # Just to prevent any collision if deBruijn graph already initialized,
        # we re-set it to empty.
        self.de_bruijn_grapth.clear()
        print('Length of the overlap graph : ',len(self.overlap_graph))
        idx_overlap_graph = 0
        print('Current date is : ',datetime.datetime.now().isoformat())
        for (key,value) in self.overlap_graph.items():
            if idx_overlap_graph == 10000 or idx_overlap_graph == 20000 or idx_overlap_graph == 30000 or idx_overlap_graph == 40000 or idx_overlap_graph == 50000 or idx_overlap_graph == 100000 :
                print('Current date is : ',datetime.datetime.now().isoformat())
                print('idx value for overlap graph is ',idx_overlap_graph,'and the length total is : ',len(self.overlap_graph))
            idx_overlap_graph+=1
            # For each key , value pair, we will build the edges of both side.
            # value can be a list of kmers.
            kmers = str(key)
            if kmers not in list_kmers_processed:
                pre_kmers = AdamskiClass.get_prefixe_kmers(kmers)
                suf_kmers = AdamskiClass.get_suffix_kmers(kmers)
                self.test_Insert_kmers_Into_DeBruijn_Graph(pre_kmers,suf_kmers)
                # Adding this kmers to the list of kmers already processed.
                list_kmers_processed.append(kmers)
            list_values = []
            list_values = list(value)
            for item in list_values:
                kmers = str(item)
                if kmers not in list_kmers_processed:
                    pre_kmers = AdamskiClass.get_prefixe_kmers(kmers)
                    suf_kmers = AdamskiClass.get_suffix_kmers(kmers)
                    self.test_Insert_kmers_Into_DeBruijn_Graph(pre_kmers,suf_kmers)
                    # Adding this kmers to the list of kmers already processed.
                    list_kmers_processed.append(kmers)

    def build_all_path(self):
        """
        To generate all the contigs find in the De Bruijn graph, we need to parse all this
        De Bruijn graph, and hence all his keys.
        """
        list_contigs    = []
        curr_str    = ''
        for key in self.de_bruijn_grapth.keys():
            for value in self.de_bruijn_grapth[key]:
                # Searching any path for each value of the current key.
                if self.de_bruijn_grapth.get(value) == None:
                    curr_str = curr_str + self.de_bruijn_grapth[key]

    def read_data_and_Build_Graph(self,nameFile,paired_reads=0,list_kmers=0):
        """
        Data from the file will look like :
             0 -> 3
             1 -> 0
             2 -> 1,6
        To insert all these data in an efficient structure, we will
        use a dictionnary.
        Adaption to the methode, to be able to read a file where we got a paired read.
        By default, it is NOT a paired_reads
        Adding a case [ with param = list_kmers ], where the file readed is filled with data such as :
        - ATG
        - ACC
        - ...
        So, by instance, for ATG , key : AT , and value : TG.
        We build the De Bruijn grapth this way.

        """
        # Clearing the DeBruijn grapth.
        self.de_bruijn_grapth.clear()
        infile = open(nameFile,'r')
        for line in infile:
            line = line.replace('\n', '')
            # Here is our adaption needed in case of paired reads.
            if paired_reads == 0 and list_kmers == 0:
                list_token = line.split(' -> ')
                key = list_token[0]
                all_values = list_token[1].split(',')
                for value in all_values:
                    self.test_Insert_kmers_Into_DeBruijn_Graph(key,value)
            elif list_kmers == 1:
                # Cases for kmers
                key = AdamskiClass.get_prefixe_kmers(line)
                value = AdamskiClass.get_suffix_kmers(line)
                self.test_Insert_kmers_Into_DeBruijn_Graph(key,value)
            else:
                list_token = line.split('|')
                # our key in this case, is the prefixe of both token.
                # As we use that as key in our dictionnary, we need an immutable structure -> tuple.
                reads_one   =   list_token[0]
                reads_two   =   list_token[1]
                key         =   AdamskiClass.get_prefixe_kmers(reads_one),AdamskiClass.get_prefixe_kmers(reads_two)
                value       =   AdamskiClass.get_suffix_kmers(reads_one),AdamskiClass.get_suffix_kmers(reads_two)
                self.test_Insert_kmers_Into_DeBruijn_Graph(key,value)

    def find_eulerian_path(self,starting_node=''):
        """
        All the edges are contained in the structure de_bruijn_graph, which is
        a dictionnary, with values representing all incident edges of the key.
        Implementation is done with link :
        www.ms.uky.ed/~lee/ma515fa10/euler.pdf
        """
        # We need a copy of de_bruijn_graph
        copy_debruijn_graph = self.de_bruijn_grapth.copy()
        # Each time, we will build a path from two adjacent edge, we will remove one value
        # of the dictionnary, so the goal would be to get an dictionnary where no anymore values
        # will remain in each key.
        eulerian_stack = []
        eulerian_path  = []
        list_incident_edges = []
        incident_edge = ''
        # By default, we can start the eulerian path by taking a randon start node
        # if no start_node specified.
        if starting_node == '':
            list_node = copy_debruijn_graph.keys()
        else:
            list_node = [starting_node]
        #for key in list(self.de_bruijn_grapth.keys()):
        for key in list_node:
            eulerian_stack.append(key)
            # The loop is only present to get one element of the dic at random.
            break
        #print self.de_bruijn_grapth.keys()
        #eulerian_stack.append(6)
        while len(eulerian_stack) > 0:
            # Searching for the eulerian path.
            top_vertex = eulerian_stack[len(eulerian_stack)-1]
            #print 'Current vertex is : ',top_vertex
            # Checking if this vertex got an incident edge not yet walked on...
            if copy_debruijn_graph.get(top_vertex) != None and len(copy_debruijn_graph[top_vertex]) > 0:
            #if len(copy_debruijn_graph[top_vertex]) > 0:
                incident_edge = list(copy_debruijn_graph[top_vertex])[0]
                #print 'All adjacent vertex of current vertex : ',copy_debruijn_graph[top_vertex]
                #print 'Adjacent vertex of current vertex is : ',incident_edge
                # Meaning there is still adjacent edge not already used.
                # there is an un-marked edge.
                # push this new edge onto the stack.
                eulerian_stack.append(incident_edge)
                # Now mark this edge -> remove this edge from the dic.
                #list_incident_edges = list(copy_debruijn_graph[top_vertex])
                #list_incident_edges.pop(0)
                edge_to_remove = list(copy_debruijn_graph[top_vertex]).pop(0)
                old_list = list(copy_debruijn_graph[top_vertex])
                old_list.pop(0)
                copy_debruijn_graph[top_vertex] = old_list
                #print 'There should be one less adjacent edges to the list : ',copy_debruijn_graph[top_vertex]
            else:
                # There is not any adjacent vertex, un-marked vertex,
                # pop the vertex off the stack, and print it to the eulerian path.
                eulerian_stack.pop(len(eulerian_stack)-1)
                eulerian_path.append(top_vertex)
        eulerian_path.reverse()
        return eulerian_path

    def print_eulerian_path(self,eulerian_list):
        """
        Printing the list the way coursera want it.
        """
        eulerian_path = ''
        for item in eulerian_list:
            eulerian_path = eulerian_path + '->' + item
        print(eulerian_path)

    def add_one_in_edges(self,list_nodes):
        """
        Adding simply 1 to the in-edges count to each elem of list_values.
        """
        # We also need to check if key exist already or not.
        for elem in list_nodes:
            if self.count_in_out_edges.get(elem) == None:
                in_edges = 1
                out_edges= 0
                self.count_in_out_edges[elem] = [in_edges,out_edges]
            else:
                in_edges = self.count_in_out_edges[elem][0] + 1
                out_edges= self.count_in_out_edges[elem][1]
                self.count_in_out_edges[elem] = [in_edges,out_edges]


    def add_in_out_edges_one_node(self,key,values):
        """
        With this element of the dictionnary, we will count his in-edges,
        and out-edges.
        """
        self.count_in_out_edges
        if self.count_in_out_edges.get(key) == None:
            # key still does not exist, creating and inserting the count value
            out_edges = len(values)
            # By default, at the beginning, we don't known in-edges,
            # so we insert simply a count of 0.
            in_edges = 0
            self.count_in_out_edges[key] = [in_edges,out_edges]
            # Updating the in-edges count for each elem in values.
            self.add_one_in_edges(values)
        else:
            # key already exist, we need to update the count value.
            in_edges = self.count_in_out_edges[key][0]
            out_edges= self.count_in_out_edges[key][1] + len(values)
            self.count_in_out_edges[key] = [in_edges,out_edges]
            # Updating the in-edges count for each elem in values.
            self.add_one_in_edges(values)

    def build_count_in_out_edge_all_nodes(self):
        """
        The function will count for each vertex in the graph, how many
        it has in-edges, and out-edges.
        This way, we will find where we have to add an edge between probably
        two nodes, and make our eulerian path, and not only the eulerian cycle.

        Again, we will use a dictionnary to count the in/out edges of a vertex.
        key/value , where value will be a list, where :
        - [0] will be the count of in-edges.
        - [1] will be the count of out-edges.
        and normally at the end, for each key in a balanced graph, [0] == [1]
        """
        for key in self.de_bruijn_grapth.keys():
            values = self.de_bruijn_grapth[key]
            self.add_in_out_edges_one_node(key,values)

    def check_unbalanced_nodes(self):
        """
        It will check for un-balanced nodes, in the De Bruijn graph.
        The structure of the De Bruijn graph, has already been builded
        and insert into : self.de_bruijn_grapth
        The checking if each nodes is balanced, is done through method
        [ self.build_count_in_out_edge_all_nodes ] , it create a new dic.
        And in this new dic, what we need to check if in-edges == out-edges
        """
        # List of nodes, where in edges are less than out edges.
        list_nodes_unbalanced_in = []
        # List of nodes, where out edges are less than out edges.
        list_nodes_unbalanced_out = []
        for key in self.count_in_out_edges.keys():
            in_edges = self.count_in_out_edges[key][0]
            out_edges= self.count_in_out_edges[key][1]
            if in_edges < out_edges:
                # This node is not balanced.
                list_nodes_unbalanced_in.append(key)
            elif in_edges > out_edges:
                list_nodes_unbalanced_out.append(key)
        return list_nodes_unbalanced_in,list_nodes_unbalanced_out


    def build_balanced_graph(self,list_nodes_unbalanced_in,list_nodes_unbalanced_out):
        """
        We got our De Bruijn graph build into dictionnary structure :
        - self.de_bruijn_grapth
        With that, we can count if each nodes is balanced [ in_edges == out_edges ],
        the function used for that is [ check_unbalanced_nodes ], it will return
        two list [ node_in , node_out ] , where node_in is the list of nodes where
        in_edges < out_edges , and the other the other way.
        So, to build a balanced graph, for each nodes in out_edges , do a link to a node
        in in_edges.
        Output  :   the last node
        """
        # So far, we do simple check, no check if one nodes got more than 2 out-going edges more
        # than in-edges.
        # We simply guess that for each out_edge there is an in_edge.
        last_node = ''
        for out_edges in list_nodes_unbalanced_out:
            for in_edges in list_nodes_unbalanced_in:
                self.test_Insert_kmers_Into_DeBruijn_Graph(out_edges,in_edges)
            last_node = out_edges
        # Returning the node who got a new out-edge. The new eulerian path, should start
        # with that one. This should be the last node of the eulerian path.
        # Don t forget that eulerian path, reverse the list as last step in the function.
        # Now, the graph present in [ self.de_bruijn_graph ] should be balanced.
        return last_node

    def string_reconstruction(self,graph_path,start_node=1,reads_paired=1,spaces_paired=0):
        """
        From the result of the eulerian path, we can re-construct the DNA genome.
        By simply, adding all of the nodes in a string object,
        first node should be full str, and the following should only be the last char of
        each node. as sub(str) == pre(str)
        Due to the fact, we also start from one special node [ start_node = 1 ], we need to
        remove the first elem of the list [ graph_path ], this elem removed is also present
        as last elem of the list.
        If start_node == 0 , no need to remove the first elem, it imply that we find only an
        eulerian cycle, and not an eulerian path. See coursera...
        """
        str_genome = ''
        # Putting the variable definition here, to be available in all function [ scope issue ]
        idx  = 0
        if start_node == 1:
            # We need to start printing the genome from the second elem [ implying removing first elem ]
            idx = 1
        else:
            idx = 0
        if reads_paired == 0:
            str_genome = str_genome + graph_path[idx]
        else:
            # Beginning printing, we get the first prefixe.
            str_genome = str_genome + graph_path[idx][0]
        idx+=1
        while idx <= len(graph_path)-1:
            # Parsing all elem of the graph path.
            # We could do all in one step, but for clarity decomposing what we are doing...
            if reads_paired == 0:
                current_elem    = graph_path[idx]
            else:
                # We first pass, we only need to get the first field of the tuple when reads_paired==1
                current_elem    = graph_path[idx][0]
            last_char       = current_elem[len(current_elem)-1]
            str_genome      = str_genome + last_char
            idx+=1
        # We reach the end of the list, and for a reads_paired, we still need to get the last suffixes elem.
        if reads_paired == 1 and spaces_paired > 0:
            # Read the forum here :
            # https://class.coursera.org/bioinformatics-001/forum/thread?thread_id=688#post-3264
            # -1 because we don't take the last elem of the list, we take it full but at the end step.
            # and still go back spaces_paried + 1 as explained in URL mentionned above.
            idx = len(graph_path) - 1 - (spaces_paired + 1)
            # idx should be sctricly lower than len because last elem shoud be taken full and not only first letter.
            while idx < len(graph_path)-1:
                # Parsing the last element, to re-construct the genome from the read pairs.
                # We now use the second elem of the tuple and not anymore the first one.
                current_elem    = graph_path[idx][1]
                first_char       = current_elem[0]
                str_genome      = str_genome + first_char
                idx+=1
        # Adding the last full elem of the list.
        str_genome  = str_genome + graph_path[idx][1]
        return str_genome

    def generate_all_binary_string(self,binary_length):
        """
        This function will generate all the binary string for a length specified
        by [ binary_length ].
        By instance :   binary_length == 2 -> 00,01,10,11
                        binary_length == 3 -> 000,001,010,100,101,011,110,111
        Output      :   a list where each elem is a binary string.

        HOWTO       :   1. the max-value of a binary string of length n is 2^n - 1
                        2. create the range of value from 0 -> 2^n - 1
                        3. convert each of these number into binary string with bin function
                        4. Cutting the string return by the bin function to not contain 0b character.
                        5. Adding the needed 0 to the string. By instance : 0001
        """
        max_bin_value   = 2**binary_length
        list_value      = range(max_bin_value)
        list_bin_str    = []
        for elem in list_value:
            # Keeping only from the second char as the first two are [ 0b ]
            cur_bin_val = bin(elem)[2:]
            if len(cur_bin_val) != binary_length:
                # We need to add the necessary 0 in front of the string.
                cur_bin_val = self.add_necessary_zeros(cur_bin_val,binary_length)
            list_bin_str.append(cur_bin_val)
        return list_bin_str

    def add_necessary_zeros(self,bin_str,bin_length):
        """
        Adding the necessary 0 in front of a binary value.
        """
        for elem in range(bin_length-len(bin_str)):
            bin_str = '0' + bin_str
        return bin_str


    @staticmethod
    def dpchange(money,coins):
        """
        Algo explained Chapter 5 , slide 71.
        coins is a list with all the piece present.
        """
        min_num_coins = numpy.zeros( (money + 1) )
        min_num_coins[0]    =   0
        for m in range(1,money+1):
            min_num_coins[m]    = sys.maxsize
            for i in range(len(coins)):
                if m >= coins[i]:
                    if min_num_coins[m-coins[i]] + 1 < min_num_coins[m]:
                        min_num_coins[m] = min_num_coins[m-coins[i]] + 1
        return min_num_coins[money]

    @staticmethod
    def manathan_tourist(n,m,down,right):
        """
        As explained in chapter 5, slide 72.

        """
        return 0

    @staticmethod
    def initialize_down_matrix_edges(namefile,row,column):
        """
        Initialize the down matrix, which means the weigth of the edges going down from nodes.
        """
        #
        matrix_edges_down   =   numpy.zeros( ( row,column+1 ) )
        infile = open(namefile,'r')
        idx_column = 0
        idx_row = 0
        for line in infile:
            line = line.replace('\n', '')
            list_token = line.split(' ')
            for weigth in list_token:
                matrix_edges_down[idx_row][idx_column] = weigth
                idx_column+=1
            idx_column = 0
            idx_row+=1
        return matrix_edges_down

    @staticmethod
    def initialize_right_matrix_edges(namefile,row,column):
        """
        Initialize the right matrix, which means the weigth of the edges going right of each nodes.
        """
        matrix_edges_right  =   numpy.zeros((row+1,column))
        infile = open(namefile,'r')
        idx_column = 0
        idx_row = 0
        for line in infile:
            print(0)

    @staticmethod
    def initialize_matrix_edges(namefile,row,column):
        matrix_edges_right  =   numpy.zeros( ( row+1,column ) )
        matrix_edges_down   =   numpy.zeros( ( row,column+1 ) )
        infile = open(namefile,'r')
        idx_column = 0
        idx_row = 0
        # We will begin processing the down matrix untill we reach the char '-' ,
        # and then processing of the right matrix will occur.
        reached_right_matrix = 0
        for line in infile:
            line = line.replace('\n', '')
            list_token = line.split(' ')
            if list_token[0] == '-':
                reached_right_matrix = 1
                idx_column = 0
                idx_row = 0
                # Continue directly reading the next line, and do not process this line.
                continue
            if reached_right_matrix == 0:
                #Processing the down matrix.
                for weigth in list_token:
                    matrix_edges_down[idx_row][idx_column] = weigth
                    idx_column+=1
                idx_column = 0
                idx_row+=1
            else:
                for weigth in list_token:
                    matrix_edges_right[idx_row][idx_column] = weigth
                    idx_column+=1
                idx_column  =   0
                idx_row+=1
        return matrix_edges_down,matrix_edges_right

    @staticmethod
    def manathan_tourist( n, m, matrix_down_edges, matrix_right_edges):
        """
        Finding the length of the longest path in this graph builded with the
        2 matrix.
        PS : range function in python output for instance :
        - range(3) = 0,1,2
        - range(1,3) = 1,2
        -> we need range(1,4) to get 1,2,3,4 .
        -> n will be n+1 in range function.
        """
        matrix_weigth = numpy.zeros( (n+1, m+1) )
        matrix_weigth[0][0] = 0
        # Processing the down edges, and computing the weigth on each node.
        for i in range(1,n+1):
            # for the matrix down edges, the value start with idx = 0 -> i-1 first reference to matrix down.
            # we use for that idx_down_matrix, to make a difference with i-1 in matrix_weigth.
            idx_down_matrix = i-1
            matrix_weigth[i][0] = matrix_weigth[i-1][0] + matrix_down_edges[idx_down_matrix][0]
        # Processing the right edges, and computing the weigth on each node.
        for j in range(1,m+1):
            # for the matrix right edges, the value start with idx = 0 -> i-1 as first ref to matrix right.
            # we use for that idx_right_matrix, to make a difference with j-1 in matrix_weigth
            idx_right_matrix = j-1
            matrix_weigth[0][j] = matrix_weigth[0][j-1] + matrix_right_edges[0][idx_right_matrix]
        for i in range(1,n+1):
            for j in range(1,m+1):
                # Again we use idx for down and right path, to make a difference with other indices.
                # This is only due to the way we build our matrix at the beginning, and the orig is 0,0 and not 1,1
                idx_down_path = i - 1
                idx_right_path= j - 1
                down_path = matrix_weigth[i-1][j] + matrix_down_edges[idx_down_path][j]
                right_path= matrix_weigth[i][j-1] + matrix_right_edges[i][idx_right_path]
                matrix_weigth[i][j] = max(down_path,right_path)
        print(matrix_weigth)
        return matrix_weigth[n][m]

    def score_edit_distance(self,str_v,str_w,indel_penalty=0):
        """
        Compute a scoring matrix needed for the edit distance.
        We need to get the minimum score [ edit distance ], between two string.
        The score is the minimum number of delete/insert/substition ine one string in order to get
        the other string.
        """
        scoring_matrix_edit_dist = numpy.zeros( ( len(str_v)+1, len(str_w)+1 ) )
        # We need to init the insert/delete score to 1 -> first column/row
        for i in range( 1,len(str_v)+1 ):
                scoring_matrix_edit_dist[i][0] = scoring_matrix_edit_dist[i-1][0] + indel_penalty
        for j in range( 1,len(str_w)+1 ):
                scoring_matrix_edit_dist[0][j] = scoring_matrix_edit_dist[0][j-1] + indel_penalty

    def get_current_score(self,list_str,idx_i,idx_j,idx_k):
        """
        In multiple longest common sub-sequence, we need to compare the character, and assign a score of 1 if all
        character are the same, and 0 otherwise.
        """
        str_v = list_str[0]
        str_w = list_str[1]
        str_x = list_str[2]
        if str_v[idx_i] == str_w[idx_j] == str_x[idx_k]:
            return 1
        else:
            return 0

    def get_max_score(self,list_str,idx_i,idx_j,idx_k):
        """
        Getting the current max score of a node, we need to search for value on the previous
        node, and the score of the current node.
        """
        str_v = list_str[0]
        str_w = list_str[1]
        str_x = list_str[2]
        # Doing the test in the same order as defined on coursera.org, varying the i , j and k idx.
        # 7 test in total
        # Backtracking direction goes to 7 direction : 1 -> 7
        current_score = -sys.maxsize
        max_score = -sys.maxsize
        current_score = self.matrix_multiple_alignment[idx_i-1][idx_j][idx_k]
        if current_score > max_score:
            max_score = current_score
            self.matrix_backtracking_multiple[idx_i][idx_j][idx_k] = 1
        current_score = self.matrix_multiple_alignment[idx_i][idx_j-1][idx_k]
        if current_score > max_score:
            max_score = current_score
            self.matrix_backtracking_multiple[idx_i][idx_j][idx_k] = 2
        current_score = self.matrix_multiple_alignment[idx_i][idx_j][idx_k-1]
        if current_score > max_score:
            max_score = current_score
            self.matrix_backtracking_multiple[idx_i][idx_j][idx_k] = 3
        current_score = self.matrix_multiple_alignment[idx_i-1][idx_j-1][idx_k]
        if current_score > max_score:
            max_score = current_score
            self.matrix_backtracking_multiple[idx_i][idx_j][idx_k] = 4
        current_score = self.matrix_multiple_alignment[idx_i-1][idx_j][idx_k-1]
        if current_score > max_score:
            max_score = current_score
            self.matrix_backtracking_multiple[idx_i][idx_j][idx_k] = 5
        current_score = self.matrix_multiple_alignment[idx_i][idx_j-1][idx_k-1]
        if current_score > max_score:
            max_score = current_score
            self.matrix_backtracking_multiple[idx_i][idx_j][idx_k] = 6
        current_score = self.matrix_multiple_alignment[idx_i-1][idx_j-1][idx_k-1] + self.get_current_score(list_str,
                                                                                                           idx_i-1,
                                                                                                           idx_j-1,
                                                                                                           idx_k-1)
        if current_score > max_score:
            max_score = current_score
            self.matrix_backtracking_multiple[idx_i][idx_j][idx_k] = 7
        return max_score



    def score_multiple_alignement(self,list_str):
        """
        Last exercice of week 7, alignement of more than 2 strings.
        """
        str_v = list_str[0]
        str_w = list_str[1]
        str_x = list_str[2]
        self.matrix_multiple_alignment = numpy.zeros( (len(str_v)+1,len(str_w)+1,len(str_x)+1 ) )
        self.matrix_backtracking_multiple = numpy.zeros( (len(str_v)+1,len(str_w)+1,len(str_x)+1 ) )
        for i in range(1, len(str_v)+1 ):
            for j in range(1, len(str_w)+1 ):
                for k in range(1, len(str_x)+1 ):
                    self.matrix_multiple_alignment[i][j][k] = self.get_max_score(list_str,i,j,k)


    def backtrack_multiple_alignment(self,list_str, idx_i, idx_j, idx_k, num_str=1):
        """
        Function to backtrack an alignement for 3 strings.
        By default, we try to backtrack printing of str 1 in list_str.
        """
        str_v = list_str[0]
        str_w = list_str[1]
        str_x = list_str[2]
        if idx_i == 0 and idx_j == 0 and idx_k == 0:
            return
        else:
            if self.matrix_backtracking_multiple[idx_i][idx_j][idx_k] == 1:
                self.backtrack_multiple_alignment(list_str,idx_i-1,idx_j,idx_k,num_str)
                print(str_v[idx_i-1])



    def lcs(self,str_v,str_w,score_matrix='',indel_penalty=0,local_alignment=False,edit_distance=False,fitting_alignment=False,overlap_alignement=False):
        """
        Implementation of algo defined in chapter 5, slide 74 from coursera.org
        We use the range function, it create our list needed.
        range(4) -> 0 1 2 3     -> we never touch an un-indexed character in the string.
        score_matrix represent the value of any match/mis-match between two amino acid.

        Edit distance Levenshtein, it is the same as global alignement but :
        - mismatch count for 1
        - match count for 0
        - insertion/deletion count for 1
        and the goal is to found the minimum at the end.
        """
        max_score_matrix = 0
        # s in the algo.
        matrix_matches_str = numpy.zeros( ( len(str_v)+1, len(str_w)+1 ) )
        # 0 represent right arrow , 1 represent diag arrow , 2 represent down arrow
        matrix_backtrack    =   numpy.zeros( ( len(str_v)+1, len(str_w)+1 ) )
        if score_matrix == '':
            # initializing the matrix with zero, to get the scoring right, as we did modity the code with matrix.
            score_matrix = numpy.zeros( ( 20, 20 ) )
        # We don't need to initialize the first row and column to zero, as it is done by numpy.zeros
        if score_matrix != '':
            # We need to init first column [ except first elem ] to 2 , down arrow
            for i in range(1,len(str_v)+1):
                matrix_backtrack[i][0] = 2
        # Only needed to put penalty on first row and column if it is a global alignment.
        if score_matrix != '' and indel_penalty != 0 and local_alignment == False:
            for i in range( 1,len(str_v)+1 ):
                matrix_matches_str[i][0] = matrix_matches_str[i-1][0] - indel_penalty
            for j in range( 1,len(str_w)+1 ):
                matrix_matches_str[0][j] = matrix_matches_str[0][j-1] - indel_penalty
        if fitting_alignment == True:
            # We need to re-set the first column to zero.
            for i in range( 1,len(str_v)+1 ):
                matrix_matches_str[i][0] = 0
        if overlap_alignement == True:
            # We need to re-set the first row and column to zero.
             for j in range( 1,len(str_w)+1 ):
                matrix_matches_str[0][j] = 0
             for i in range( 1,len(str_v)+1 ):
                matrix_matches_str[i][0] = 0
        for i in range(1, len(str_v)+1 ):
            for j in range(1, len(str_w)+1 ):
                # Init these two var here, to be available somewhere below...
                score_match = 0
                score_missmatch = 0
                previous_i = matrix_matches_str[i-1][j] - indel_penalty
                previous_j = matrix_matches_str[i][j-1] - indel_penalty
                # match between two amino acid in the string v and w.
                if str_v[i-1] == str_w[j-1]:
                    # Getting the idx in the matrix with the amino acid.
                    idx_matrix = self.idx_matrix_amino_acid[str_v[i-1]]
                    # Getting the score of this match in the matrix.
                    score_match = score_matrix[idx_matrix][idx_matrix]
                    previous_i_j_diag = matrix_matches_str[i-1][j-1] + score_match
                    if local_alignment == True:
                        if edit_distance == False:
                            matrix_matches_str[i][j] = max(0,previous_i,previous_j,previous_i_j_diag)
                        else:
                            matrix_matches_str[i][j] = min(previous_i,previous_j,previous_i_j_diag)
                        # We need to store where is the max value in the matrix.
                        idx_i = self.max_value_matrix_score[0]
                        idx_j = self.max_value_matrix_score[1]
                        max_score_matrix = matrix_matches_str[idx_i][idx_j]
                        if matrix_matches_str[i][j] >= max_score_matrix:
                            self.max_value_matrix_score = (i,j)
                    else:
                        if edit_distance == False:
                            matrix_matches_str[i][j] = max(previous_i,previous_j,previous_i_j_diag)
                        else:
                            matrix_matches_str[i][j] = min(previous_i,previous_j,previous_i_j_diag)
                else:
                    # Getting the idx in the matrix with the amino acid.
                    idx_matrix_v = self.idx_matrix_amino_acid[str_v[i-1]]
                    idx_matrix_w = self.idx_matrix_amino_acid[str_w[j-1]]
                    # Getting the score of this match in the matrix.
                    score_missmatch = score_matrix[idx_matrix_v][idx_matrix_w]
                    previous_i_j_diag = matrix_matches_str[i-1][j-1] + score_missmatch
                    if local_alignment == True:
                        if edit_distance == False:
                            matrix_matches_str[i][j]    =   max(0,previous_i,previous_j,previous_i_j_diag)
                        else:
                            matrix_matches_str[i][j] = min(previous_i,previous_j,previous_i_j_diag)
                        # We need to store where is the max value in the matrix.
                        idx_i = self.max_value_matrix_score[0]
                        idx_j = self.max_value_matrix_score[1]
                        max_score_matrix = matrix_matches_str[idx_i][idx_j]
                        if matrix_matches_str[i][j] >= max_score_matrix:
                            self.max_value_matrix_score = (i,j)
                    else:
                        if edit_distance == False:
                            matrix_matches_str[i][j]    =   max(previous_i,previous_j,previous_i_j_diag)
                        else:
                            matrix_matches_str[i][j] = min(previous_i,previous_j,previous_i_j_diag)
                if indel_penalty == 0:
                    if matrix_matches_str[i][j] == matrix_matches_str[i-1][j]:
                        matrix_backtrack[i][j] = 2
                    elif matrix_matches_str[i][j] == matrix_matches_str[i][j-1]:
                        matrix_backtrack[i][j] = 0
                    elif matrix_matches_str[i][j] == ( matrix_matches_str[i-1][j-1] + 1 ):
                        matrix_backtrack[i][j] = 1
                else:
                    # Using the code with initialized matrix and indel_penalty set.
                    if matrix_matches_str[i][j] == matrix_matches_str[i-1][j] - indel_penalty:
                        # Best score is to get down
                        matrix_backtrack[i][j] = 2
                    elif matrix_matches_str[i][j] == matrix_matches_str[i][j-1] - indel_penalty:
                        matrix_backtrack[i][j] = 0
                    elif matrix_matches_str[i][j] == ( matrix_matches_str[i-1][j-1] + score_match ) or matrix_matches_str[i][j] == ( matrix_matches_str[i-1][j-1] + score_missmatch ):
                        # We should print the char, even if it does not match, it is a better score than right arrow
                        # or down arrow.
                        matrix_backtrack[i][j] = 1
                    elif matrix_matches_str[i][j] == 0:
                        # The best predecessor is the source node.
                        matrix_backtrack[i][j] = 3
        # Returning the last value created -> len(str) - 1
        #return matrix_matches_str[len(str_v)-1][len(str_w)-1],matrix_backtrack
        return matrix_matches_str,matrix_backtrack


    def output_lcs(self,matrix_backtrack,str_v,i,j,global_alignement=False,first_str=True,fitting_alignment=False):
        """
        PS : We should check if the function still works properly for a normal alignement.
                Because this one work now for global alignement.

        As said in our previous function lcs, in backtrack matrix :
        - 0 represent , right arrow
        - 1 represent , diag arrow
        - 2 represent , down arrow

        PS : As normally, we do the alignement against two string, the argument first_str means that we need to print
                this alignement against the first str, otherwise the second.
        """
        if fitting_alignment == False:
            if i == 0 and j == 0:
                return
        else:
            # Depending which string is the shorted, we will stop when we encounter i=0 or j=0.
            if i == 0 or j == 0:
                return
        if matrix_backtrack[i][j] == 2:
            self.output_lcs(matrix_backtrack,str_v,i - 1,j,global_alignement,first_str,fitting_alignment)
            if global_alignement == True:
                if first_str == True:
                    print(str_v[i-1],)
                else:
                    print('-',)
        elif matrix_backtrack[i][j] == 0:
            self.output_lcs(matrix_backtrack,str_v,i,j-1,global_alignement,first_str,fitting_alignment)
            if global_alignement == True:
                if first_str == True:
                    print('-',)
                else:
                    print(str_v[j-1],)
        else:
            self.output_lcs(matrix_backtrack,str_v,i-1,j-1,global_alignement,first_str,fitting_alignment)
            # Start printing from i-1 instead of i, because the string start
            # at position 1, and the matrix start at position 0
            if global_alignement == True:
                if first_str == True:
                    print(str_v[i-1],)
                else:
                    print(str_v[j-1],)
            else:
                print(str_v[j-1],)

    def test_and_insert_predecessor_node(self,node,predecessor,weight):
        if self.predecessor_nodes.get(node) == None:
            # This is a new node.
            self.predecessor_nodes[node] = [(predecessor,weight)]
        else:
            # node already exist, inserting a nez value into the list.
            cur_list = self.predecessor_nodes[node]
            cur_list.append((predecessor,weight))
            self.predecessor_nodes[node] = cur_list

    def read_data_weigth_edges(self,namefile):
        """
        It will read the data in a specific file. The data in the file are :
        0->1:7
        0->2:4
        2->3:2
        """
        self.predecessor_nodes.clear()
        infile = open(namefile,'r')
        for line in infile:
            line = line.replace('\n', '')
            list_token = line.split('->')
            node_orig = list_token[0]
            destination = list_token[1].split(':')
            node_dest = destination[0]
            weigth_to_dest = destination[1]
            if self.edges_weigth_dag.get(node_orig) == None:
                # This is a new node -> creating it...
                self.edges_weigth_dag[node_orig] = [(node_dest,weigth_to_dest)]
            else:
                # this node already exist, adding a new tuple to the list of destination.
                list_dest = self.edges_weigth_dag[node_orig]
                list_dest.append( (node_dest,weigth_to_dest))
                self.edges_weigth_dag[node_orig] = list_dest
            # Creating the dictionnary for the relation of predecessor nodes.
            self.test_and_insert_predecessor_node(node_dest,node_orig,weigth_to_dest)

    def get_weight_nodes(self,node,start_node):
        """
        Function name say it all.
        We will stop if we encounter the start_node, or if we goes through this node, and it is not the
        biggest value, we will force to take it.
        """
        weight_node = 0
        max_node_weight = 0
        path_max_node = ''
        if self.weight_nodes.get(node) == None:
            # Still not yet a weight of node -> finding it recursively
            #print 'Getting node : ',node
            #print '     The predecessor are : ',self.predecessor_nodes[node]
            for tuple_predecessor in self.predecessor_nodes[node]:
                node_name = tuple_predecessor[0]
                weight_edge = tuple_predecessor[1]
                node_value_prede = self.get_weight_nodes(node_name,start_node)
                weight_node = int(node_value_prede) + int(weight_edge)
                if weight_node > max_node_weight:
                    max_node_weight = weight_node
                    path_max_node = node_name
                if node_name == start_node:
                    max_node_weight = weight_node
                    path_max_node = node_name
                    break
            self.weight_nodes[node] = (max_node_weight,path_max_node)
            self.node_used.append(path_max_node)
            print('For node : ',node,' The max value is : ', max_node_weight ,'which goes throught node : ',path_max_node)
            return max_node_weight
        else:
            return self.weight_nodes[node][0]

    def init_start_node(self):
        for node in self.edges_weigth_dag.keys():
            if self.predecessor_nodes.get(node) == None:
                # We got a starting node in the graph, no predecessor -> weight node = 0
                self.weight_nodes[node] = (0,None)

    def build_dag_grath(self,start_node,sink_node):
        """
        Only used for the exercice Chapter 5 , slide 74.

        We need to loop over all the key we inserted, and from that we build the graph.
        Each node in the graph, got :
        - the destination where he can go
        - the node it came from
        - his current weigth [ max weigth of the edge of all of his predecessor ]
        """
        # Clearing the graph just for safery.
        self.graph_dag.clear()
        list_dest_node = self.edges_weigth_dag[start_node]
        weight_start_node = 0
        max_weight = 0
        orig_node = ''
        for dest in list_dest_node:
            dest_node = dest[0]
            weight_node = dest[1]
            if ( weight_node + weight_start_node ) > max_weight:
                max_weight = weight_node + weight_start_node
                orig_node = start_node


        list_nodes = []
        list_nodes = self.edges_weigth_dag.keys()
        while len(list_nodes) != 0:
            # We still have some nodes to processes.
            return 1

    def build_scoring_matrix(self,namefile):
        """
        Chapter 5 , slide 74. We need a matrix to score each LCS.
        The matrix is in a file
        """
        # Initializing a matrix with the score.
        matrix_score = numpy.zeros( (20,20) )
        infile = open(namefile, 'r')
        listToken = []
        # We will fill the matrix, row by row with this list.
        value_to_insert = []
        idx_amino_acid = 0
        idx_column_matrix = 0
        for line in infile:
            line = line.replace('\n', '')
            listToken = line.split()
            # The first line should be the amino acid, we remove it.
            # Each time, the first char is the amino acid.
            amino_acid = listToken[0]
            self.idx_matrix_amino_acid[amino_acid] = idx_amino_acid
            # Because they may be more than 1 space between values, we use the simple split(), which handle that.

            matrix_score[idx_amino_acid] = listToken[1:]
            # Parsing the rest of the line, and getting the scoring against two amino acid
            """
            for score in all_values:
                print score,'-next-'
                matrix_score[idx_amino_acid] = all_values
                idx_column_matrix+=1
            """
            # Resetting the idx for column, for the next line.
            idx_column_matrix = 0
            # Incrementing the row index for the matrix.
            idx_amino_acid+=1
        return matrix_score

    def backtrack_local_alignement(self,matrix_score,idx_i,idx_j,str_v,str_w,indel_penalty,matrix_score_matches,first_str=True):
        """
        In matrix_score, there will be all the score build with the function lcs.
        And in a classe variable [ max_value_matrix_score ] , will contains the idx i and j
        of the max value of the matrix. This is needed, because from there, we will begin
        printing the local alignment as asked in coursera.org
        We just have to pay attention, on the way we put the score in the matrix build in function lcs.
        If it start from righ , diag , down, or any other order.
        """
        if matrix_score[idx_i][idx_j] == 0:
            return
        else:
            # We first need to get the score of match/missmatch of two character in the string v and w, to known
            # where we did go
            idx_matrix_v = self.idx_matrix_amino_acid[str_v[idx_i-1]]
            idx_matrix_w = self.idx_matrix_amino_acid[str_w[idx_j-1]]
            # Getting the score of this match in the matrix.
            score = matrix_score_matches[idx_matrix_v][idx_matrix_w]
            previous_i_j_diag = matrix_score[idx_i-1][idx_j-1] + score
            if matrix_score[idx_i][idx_j] == matrix_score[idx_i-1][idx_j] - indel_penalty:
                self.backtrack_local_alignement(matrix_score,idx_i-1,idx_j,str_v,str_w,indel_penalty,matrix_score_matches,first_str)
                if first_str == True:
                    print(str_v[idx_i-1],)
                else:
                    print('-',)
            elif matrix_score[idx_i][idx_j] == matrix_score[idx_i][idx_j-1] - indel_penalty:
                self.backtrack_local_alignement(matrix_score,idx_i,idx_j-1,str_v,str_w,indel_penalty,matrix_score_matches,first_str)
                if first_str == True:
                    print('-',)
                else:
                    print(str_w[idx_j-1],)
            elif matrix_score[idx_i][idx_j] == previous_i_j_diag:
                self.backtrack_local_alignement(matrix_score,idx_i-1,idx_j-1,str_v,str_w,indel_penalty,matrix_score_matches,first_str)
                if str_v[idx_i-1] == str_w[idx_j-1]:
                    # This is a match, we can print whatever character of both string, str_v or str_w
                    print(str_v[idx_i-1],)
                else:
                    if first_str == True:
                        print(str_v[idx_i-1],)
                    else:
                        print(str_w[idx_j-1],)

    @staticmethod
    def get_idx_max_value_last_column_matrix(matrix):
        """
        This method is only for execice fitting alignement.
        To print the string, we need to get the max value from the last column of the matrix.
        And then from there, do the global alignment.
        """
        # Getting the row idx from each column of the matrix, where there is the max value.
        list_row_idx = matrix.argmax(axis=0)
        # As the matrix should be a 2 dimensions. 1 should reference the number of column in matrix
        num_column = matrix.shape[1]
        return (list_row_idx[num_column-1],num_column-1)

    @staticmethod
    def get_idx_max_value_last_row_matrix(matrix):
        """
        This method is only for exercice overlap alignment. Because we need to known where is
        the max value on the last row of the matrix. And on this last row, we need to get the last value found
        , and not the first found on this row.
        """
        num_row = matrix.shape[0]
        # Getting all elem of the last row, except the zero put in the first column.
        all_elem = matrix[num_row-1][1:]
        # As by default, all max function seems to found the first max value, we will implement our max value
        # because in our case, we need the last max value.
        idx_j = 1
        max_value = -sys.maxsize
        idx_j_max_value = 1
        for elem in all_elem:
            if elem >= max_value:
                max_value = elem
                idx_j_max_value = idx_j
            idx_j+=1
        return (num_row-1,idx_j_max_value)

    @staticmethod
    def get_middle_edge(matrix_score,len_str_w):
        """
        Finding the row and column number where the middle edge is.
        """
        middle = len_str_w/2
        # Getting the row idx from each column of the matrix, where there is the max value.
        list_row_idx = matrix_score.argmax(axis=0)
        return (list_row_idx[middle],middle)

    def build_back_list(self,list_orig,list_inserting,begin_idx,end_idx):
        """
        Re-building the original list.
        There is 3 case possible, idx_x = 0 , idx_x > 0 and idx_x < len(list) , otherwise idx_x == len(list)-1
        it is the latest item in the list.
        """
        if begin_idx == 0:
            return list_inserting+list_orig[end_idx:]
        elif end_idx < len(list_orig):
            # need to rebuild, before and after the new list.
            return list_orig[0:begin_idx]+list_inserting+list_orig[end_idx:]
        else:
            # We reach the last item of the list, adding the new list at the end.
            return list_orig[0:begin_idx]+list_inserting


    def greeding_sorting(self,list):
        """
        Greedy sorting algo as described in coursera.org for bio-informatics.
        """
        approx_reversal_distance = 0
        cur_idx = 0
        rebuilded_list = list[:]
        for x in range(1,len(list)+1):
            # We do not known if it is x or -x which is present in the list -> need to check.
            if rebuilded_list.count(x) > 0:
                num_to_handle = x
            elif rebuilded_list.count(-x) > 0:
                num_to_handle = -x
            else:
                print("Number not found in the list, not GOOD")
                print('The number is : ',x)
                sys.exit()
            # doing processing on a positive number
            # We can make an improvement, by searching only from x -> till end of list, and not always all the list.
            idx_x = rebuilded_list.index(num_to_handle)
            if cur_idx != idx_x:
                # We need to sort, as the element is not on the right place.
                # reversing from the current idx -> the number found in the str [ +1 which take the number it-self also ]
                list_to_reverse = rebuilded_list[cur_idx:idx_x+1]
                new_list_builded = self.sorting_reversal(list_to_reverse)
                # Re-build the original list with this new list sorted_reversal.
                # need to do some re-building...
                # There is 3 case possible, idx_x = 0 , idx_x > 0 and idx_x < len(list) , otherwise idx_x == len(list)-1
                # it is the latest item in the list.
                rebuilded_list = self.build_back_list(rebuilded_list,new_list_builded,cur_idx,idx_x+1)
                self.printing_list(rebuilded_list)
                approx_reversal_distance+=1
            if rebuilded_list[x-1] < 0:
                rebuilded_list[x-1] = -rebuilded_list[x-1]
                self.printing_list(rebuilded_list)
                approx_reversal_distance+=1
            cur_idx+=1
            # normally, we shouldn't use this variable cur_idx, as it can be deducted from x instead.
        return approx_reversal_distance

    def printing_list(self,list):
        """
        Needed to output the list in the proper format by coursera, meaning a + in front of number if positif.
        """
        str_to_print = ''
        str_to_print = str_to_print+'('
        for x in list:
            if x > 0:
                str_to_print = str_to_print+'+'+str(x)+' '
            else:
                str_to_print = str_to_print+str(x)+' '
        # Removing the last whitespace inserted.
        str_whithout_last_space = str_to_print.rstrip()
        str_whithout_last_space = str_whithout_last_space+')'
        print(str_whithout_last_space)

    def sorting_reversal(self,list_num):
        """
        Function used to reverse a list, and inverting sign of each number present in the list.
        """
        # Reversing the list
        list_reversed = list_num[::-1]
        # Building a matrix of 1 dim.
        matrix = numpy.array( list_reversed )
        # Inverting sign of all number in the list
        return (matrix*-1).tolist()

    @staticmethod
    def compute_GC_content( dna):
        """
            As defined on rosalind.org , the GC content of a DNA, is the ratio of C or G over
            the whole dna string.
            GC Ratio for :
            ATCG        -> 50.0
            ACG         -> 66.6
            AATTCG      -> 33.3
        """
        return ((dna.count('C')+dna.count('G'))/float(len(dna)))*100


    @staticmethod
    def load_fasta_files(input_file, unix_file=True):
        """
            The function will return a dictionnary containing all the data in the file.
            The key will be the FASTA ID, and the content his DNA associated.
            Example :
                >Rosalind_6404
                CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
                TCCCACTAATAATTCTGAGG

                -> key          : Rosalind_6404
                   key content  :    CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
                                     TCCCACTAATAATTCTGAGG
        """
        last_id_fasta = ''
        id_fasta = ''
        dic_fasta = {}
        current_dna = ''
        infile = open( input_file,'r')
        for line in infile:
            if unix_file:
                line = line.replace('\n', '')
            else:
                print('### TO DO ###')
                print('Need to check what is the char end of line for another OS...')
            if line[0] == '>':
                id_fasta = line[1:]
                if current_dna != '':
                    # Inserting a new key + value into the dictionnary.
                    # We don t use id_fasta because it just has been overwriten with
                    # current reading line. That s why keeping the last reference.
                    dic_fasta[last_id_fasta] = current_dna
                    # Re-initializing the dna for the next entry.
                    current_dna = ''
            else:
                last_id_fasta = id_fasta
                # This is the dna string, continuing building it...
                current_dna = current_dna + line
        # Inserting the last dna and fasta id into the dic.
        dic_fasta[id_fasta] = current_dna
        return dic_fasta

    def compute_breakpoint_number(self,list_permutation):
        """
        Finding the number of breakpoint present in the list.
        """
        number_break = 0
        adjacency = 0
        if list_permutation[0] != 1:
            # Because we have implicitely 0 in front of the list
            number_break+=1
        if list_permutation[len(list_permutation)-1] != len(list_permutation):
            # Because we have implicitely the last number at the end of the list.
            number_break+=1
        idx = 0
        while idx < len(list_permutation)-2:
            if list_permutation[idx]+1 == list_permutation[idx+1]:
                adjacency+=1
            else:
                number_break+=1
            idx+=1
        return number_break

    def build_tries_construction(self, list_pattern):
        """
        Building the list of adjancy node for the tries contruction.
        """
        # We first need to initialize the tries construction, the root node with 1.
        first_pattern = list_pattern[0]
        rest_list_pattern = list_pattern[1:]
        first_str_pattern = first_pattern[0]
        rest_first_str_pattern = first_pattern[1:]
        next_idx = self.index_tries_construction + 1
        self.index_tries_construction = next_idx
        self.tries_construction[1] = [(next_idx,first_str_pattern)]
        # We just build the first elem in the tries construction.
        # And we now need to continue the construction of this first pattern.
        for x in rest_first_str_pattern:
            next_idx = self.index_tries_construction + 1
            self.tries_construction[self.index_tries_construction] = [(next_idx,x)]
            self.index_tries_construction = next_idx
        # We now need to construct the rest of the tree (tries) with all the other pattern in list_pattern.
        for pattern in rest_list_pattern:
            # we need to found where to insert in the tries.
            # We always start to insert a pattern from the root node -> node number 1
            idx = self.found_idx_in_tries(pattern,1)
            # How deep did we go in the tries ? This will tell us how much of the rest of the pattern
            # need to be inserted in the tries.
            deep_level = self.deep_level_in_tries
            rest_pattern_to_be_inserted = pattern[deep_level:]
            self.insert_pattern_in_tries(rest_pattern_to_be_inserted,idx)
            # We need to reset the level of deep we goes through last time.
            self.deep_level_in_tries = 0


    def insert_pattern_in_tries(self, pattern, node_number):
        """
        The function will insert the rest of pattern in the tries, the beginning of the
        pattern match may be a branch in the tries, that's why we only insert a part of the pattern.
        Insertion will begin from node number [ node_number ]
        """
        # First checking if the node to insert is a new one or not ?
        # Getting first the list of letter already in the tries.
        if self.tries_construction.get(node_number) == None:
            # this node number does not exist yet... creating...
            next_idx = self.index_tries_construction + 1
            self.index_tries_construction = next_idx
            letter = pattern[0]
            self.tries_construction[node_number] = [(next_idx,letter)]
            # Checking if there is still something to be inserted...
            if len(pattern[1:]) > 0:
                # still something to be inserted in the tries...
                self.insert_pattern_in_tries(pattern[1:],next_idx)
        else:
            list_letter = self.tries_construction[node_number]
            next_idx = self.index_tries_construction + 1
            self.index_tries_construction = next_idx
            letter = pattern[0]
            list_letter.append( (next_idx,letter) )
            # Inserting this new branch to the list of the node.
            self.tries_construction[node_number] = list_letter
            if len(pattern[1:]) > 0:
                # still something to be inserted in the tries...
                self.insert_pattern_in_tries(pattern[1:],next_idx)

    def found_idx_in_tries(self, pattern, node_number):
        """
        The first time we call this function, it should be with node_number = 1 , 1 is the root node.
        Searching recursively in all node of the tries, where to insert the current letter.
        We now got a tries [ tree ] , and we need to traverse the tree to found
        where we can start inserting a new pattern.
        If the letter is not yet present in the tries, we insert from root node [ 1 ], otherwise
        somewhere below in the tries.
        """
        list_from_cur_node = []
        list_from_cur_node = self.tries_construction[node_number]
        # node number does not exist, and won't never exist.
        found_from_node_number = node_number
        for elem in list_from_cur_node:
            if elem[1] == pattern[0]:
                # Continue deper in the tries and searching for continuing pattern matching.
                found_from_node_number = self.found_idx_in_tries(pattern[1:],elem[0])
                self.deep_level_in_tries = self.deep_level_in_tries + 1
                # We break the for loop, because ASAP that we found a match, we go deper in the tries...
                break
        return found_from_node_number

    def print_tries_structure(self):
        """
        This function will print the tries structure as the way bio-informatica coursera.org want it, to
        get a good result, and hence the point associated with it.
        """
        for key in self.tries_construction.keys():
            # Iterating over all node of the tries structure.
            for elem in self.tries_construction[key]:
                # Iterating over each letter present in the current node of the tries.
                print(key, elem[0], elem[1])

    def prefix_trie_matching(self, text, idx_match_search):
        """
        Algo as defined by coursera.org with the trie structure.
        """
        cur_idx = 0
        symbol = text[cur_idx]
        cur_idx+=1
        v = 1
        pattern_spelled = ''
        while True:
            #print 'parameter passed : ',v,symbol
            edge_found,node_number = self.is_edge_symbol_found(v,symbol)
            #print edge_found,node_number
            if edge_found == True:
                pattern_spelled = pattern_spelled + symbol
                if len(text[cur_idx:]) > 0:
                    symbol = text[cur_idx]
                else:
                    # as end of str, just to init symbol to something...
                    symbol = ''
                cur_idx+=1
                v = node_number
            elif self.tries_construction.get(v) == None:
                # print pattern_spelled, idx_match_search,
                # We only need the idx where the pattern does exist in the string...
                print(idx_match_search,)
                return
            else:
                # print 'Not match found'
                # We do not need to print this, output will be easier for coursera.org
                return

    def trie_matching(self, text):
        """
        Algo as defined on coursera.org for trie structure.
        """
        idx = 0
        while len(text[idx:]) > 0:
            self.prefix_trie_matching(text[idx:],idx)
            idx+=1


    def is_edge_symbol_found(self, node_number, letter):
        """
        The function will only tell us if at node [ node_number ] there is a letter [ letter ]
         meaning an edge v -> w with letter [ letter ]
        """
        if self.tries_construction.get(node_number) == None:
            return False,False
        else:
            all_edges = self.tries_construction[node_number]
            for edge in all_edges:
                if edge[1] == letter:
                    return True,edge[0]
                    # Do not need to process further the list of edges...
                    break
        return False,False

    def get_all_suffix(self, text):
        """
        Will build a list with all suffixe of text, including the '$' character.
        """
        text = text + '$'
        list_suffixe = []
        idx = 0
        while idx < len(text):
            list_suffixe.append(text[idx:])
            idx+=1
        return list_suffixe


    def search_longest_repeat(self):
        """
        some text to add ...
        """
        list_root_nodes = self.tries_construction[1]
        longest_pattern = ''
        for node in list_root_nodes:
            self.search_recur_longest_repeat(node, node[1], longest_pattern)


    def search_recur_longest_repeat(self, cur_node, letter, longest_pattern):
        """
        some text to add...
        """
        node_reference = cur_node[0]
        if self.tries_construction.get(node_reference) == None:
            # We did reach a leaf, we do not processes further...
            return
        list_sub_nodes = self.tries_construction[node_reference]
        new_longest_pattern = longest_pattern + letter
        if len(list_sub_nodes) > 1:
            # We got multiple sub-node, and hence a potential longest repeat sequence.
            if len(new_longest_pattern) > len(self.longest_repeat):
                self.longest_repeat = new_longest_pattern
        for node in list_sub_nodes:
            self.search_recur_longest_repeat(node, node[1], new_longest_pattern)


    def merge_non_branching_path(self, current_node, reference_node_number, pattern):
        """
        NOT YET USED AS FUNCTION
        This function will build up the suffix tree from the tries structure.
        So far, it won't be very fast... need to improve that !
        """
        if self.tries_construction.get(current_node) == None:
            # We got a leaf, nothing more to be done...
            # Re-building the potential non branching tree...
            self.tries_construction[reference_node_number] = [(current_node, pattern)]
        else:
            list_sub_nodes = self.tries_construction[current_node]
            if len(list_sub_nodes) == 1:
                # We need to go further down the tree, and continue re-assembling the pattern...
                # As list_sub_node should be a list with only one tuple.
                node = list_sub_nodes[0]
                pattern = pattern + node[1]
                self.merge_non_branching_path(node[0], reference_node_number, pattern)
            elif len(list_sub_nodes) > 1:
                # Rebuilding a non-branching path.
                print(0)
            else:
                print("We should never PRINT THIS LINE !!!!!!!!!!!")

    def get_longest_prefixe_match(self, pattern1, pattern2):
        """
        We return the longest commumn prefixe betwwen the two patter [ str ]
        pattern2 should be the longest pattern.
        pattern1 the one already in the tree.
        """
        idx = 0
        commun_pattern = ''
        while idx < len(pattern2) and idx < len(pattern1):
            if char1[idx] == char2[idx]:
                idx+=1
                commun_pattern = commun_pattern + char1
            else:
                break
        return commun_pattern

    def get_cyclic_rotation(self, text):
        """
        The function will output a list containing all the strings of text,
        which represent a cyclic rotation of the original text.
        """
        list_cyclic = []
        list_cyclic.append(text)
        idx = len(text)-1
        rotated_text = ''
        while ( idx > 0 ):
            rotated_text = text[idx:] + text[0:idx]
            list_cyclic.append(rotated_text)
            idx = idx - 1
        return list_cyclic

    def get_BWT( self, list_text):
        """
        Implement the Burrows-Wheeler transform, given a matrix with all lines
        of text.
        """
        # the sorting occur in place -> no new for a new variable.
        list_text.sort()
        # All items in the list will have the same length.
        idx_last_item = len(list_text[0])-1
        bw_transform = ''
        for item in list_text:
            bw_transform = bw_transform + item[idx_last_item]
        return bw_transform

    def get_BW_transform( self, text):
        """
        Easier for us, if all the Burrows-Wheeler occurs in only one function.
        """
        list_cyclic_rotation_text = self.get_cyclic_rotation(text)
        return self.get_BWT( list_cyclic_rotation_text )

    def fill_matrix_BWT(self , text, fill_first_row=True, do_preprocessing=True):
        """
        With the BW transform of a text, we can fill a matrix, and find back the original
        text, as described on coursera.org
        """
        # As we cannot use numpy with character item in it.
        # We will use a list of list
        matrix_BW = range(len(text))
        idx = 0
        while idx < len(text):
            matrix_BW[idx] = range(len(text))
            idx+=1
        # Initializing the last column of the matrix with the BWT.
        idx_y = 0
        idx_x = len(text)-1
        for char in text:
            matrix_BW[idx_x][idx_y] = char
            idx_y+=1
        # Init the first column of the matrix with text re-ordered lexicographically
        orig_text = text[:]
        sorted_text = self.transform_text_to_list_sorted( text)
        idx_x = 0
        idx_y = 0
        for char in sorted_text:
            matrix_BW[idx_x][idx_y] = char
            idx_y+=1
        if do_preprocessing == True:
            # Pre-process the matrix, to known in advance when we got the character from the first column,
            # if it is the first, second, or third, ... a char in this list.
            self.pre_process_first_column_BW_matrix(matrix_BW)
            # Now the real re-construction of the full matrix will begin.
            # and finally getting the original string.
            # Adaption has been done on get_x_elem, to use now a dictionnary, instead of always going through
            # the list each time.
            self.pre_process_last_column_BW_matrix( matrix_BW)
        if fill_first_row == True:
            self.fill_first_row_matrix(matrix_BW)
        return matrix_BW

    def get_inverse_BW_transform(self, text):
        """
            As input we get the string, which< is the BW transform of orig_text, and we need to
            find this orig_text which did generate the text.
        """
        # Fill the first and last column of the matrix, and also the first row of the matrix, which is
        # our original text with $ in front.
        matrix_BW = self.fill_matrix_BWT( text)
        # Re-building the string from the first row.
        orig_text = ''
        # starting from 1, this way the '$' is not in front of the text, we will add it at the end.
        idx = 1
        while idx < len(matrix_BW[0]):
            orig_text = orig_text + matrix_BW[idx][0]
            idx+=1
        orig_text = orig_text + '$'
        return orig_text

    def fill_first_row_matrix(self, matrix_BW):
        """
            Now that the last and first column of the matrix BW is filled, we should be able
            to re-construct the first row, which is actually, the string that did make the matrix.
        """
        idx_last_column_BW = len(matrix_BW[0])-1
        index_char = matrix_BW[idx_last_column_BW].index('$')
        # tup will represent the tuple which are present in the first column of the matrix.
        tup = matrix_BW[0][index_char]
        char_first_column = tup[0]
        # If this char is not the first on the column...
        count_char = tup[1]
        # Re-setting the first char of the matrix, to a single char and not a tuple, it will be
        # easier afterwards, to re-compose the string which made the matrix.
        matrix_BW[0][0] = '$'
        matrix_BW[1][0] = char_first_column
        # matrix_BW is a square matrix, we could have choose len(matrix_BW[1])...
        for idx_row in range(2, len(matrix_BW[0])):
            index_char = self.get_x_elem(matrix_BW[idx_last_column_BW], char_first_column, count_char)
            tup = matrix_BW[0][index_char]
            char_first_column = tup[0]
            count_char = tup[1]
            matrix_BW[idx_row][0] = char_first_column

    def pre_process_last_column_BW_matrix(self, matrix_BW):
        """
            As the function get_x_elem is executed each times we need an elem from the last column and moreover,
            as there is more than 10.000 element. The total run time is way too much...
            So, we will build a dictionnary containing the index where a character and his count is located on the last column.
            By instance :
                dic[(a, 3)] = 4
                meaning that the third 'a' character in the last column of the matrix is located at index 4.
        """
        count_A = 1
        count_T = 1
        count_C = 1
        count_G = 1
        idx = 0
        # strange way to get the last idx, but it is working...
        idx_last_column = len(matrix_BW[0])-1
        for char in matrix_BW[idx_last_column]:
            if self.idx_character_last_column.get((char, 1)) == None:
                self.idx_character_last_column[(char, 1)] = idx
            else:
                if char == 'A':
                    count_A+=1
                    self.idx_character_last_column[(char, count_A)] = idx
                elif char == 'T':
                    count_T+=1
                    self.idx_character_last_column[(char, count_T)] = idx
                elif char == 'C':
                    count_C+=1
                    self.idx_character_last_column[(char, count_C)] = idx
                elif char == 'G':
                    count_G+=1
                    self.idx_character_last_column[(char, count_G)] = idx
            idx+=1

    def pre_process_first_column_BW_matrix(self, matrix_BW):
        """
        Instead of just placing the character in the first column of the matrix BW, we will put a tuple,
        which contains, the character it-self, and if it is the first, second, third, ... character of this list.
        It makes our life easier, when we have to find the match on the last column of matrix BW.
        """
        dic_count_character = {}
        idx = 0
        for char in matrix_BW[0]:
            if dic_count_character.get(char) == None:
                # We got a new character, inserting in the dictionnary, and replacing this item in the matrix.
                dic_count_character[char] = 1
                matrix_BW[0][idx] = (char, 1)
            else:
                # There was already one of the current char in the dic
                count_char = dic_count_character[char]
                count_char+=1
                dic_count_character[char] = count_char
                matrix_BW[0][idx] = (char, count_char)
            idx+=1

    def get_x_elem(self, list, character, number):
        """
            We now use a dictionnary, to known where the Xth character is located on the last column.
            Previously, we were searching through the list [ list ] , but as there was more than 10.000 in the list,
            and we needed to do a lot of search, it took too much time...
        """
        return self.idx_character_last_column[(character, number)]

    def transform_text_to_list_sorted(self, text):
        """
            As we cannot sort directly a string, we will first put all the string
            in a list, and then sort this list.
        """
        list = []
        for x in text:
            list.append(x)
        list.sort()
        return list

    def build_suffixe_array(self, text, partial=False, k=0):
        """
            Building a suffixe array.
            We will implement it with a list, where each elem if the list, is a tuple which containt
                1. the suffixe it-self
                2. the index where this suffixe is located in the text.
        """
        # TO DO build the suffixe list inside the loop while , because otherwise take too much memory.
        idx = 0
        suffixe_array = []
        suffixe = ''
        while idx < len(text):
            # using buffer method is way much faster on big string.
            suffixe = buffer(text, idx, len(text))
            #suffixe = text[idx:]
            # Creating a new tuple with suffixe, and position where this suffixe is located in original text.
            new_tup = (suffixe, idx)
            suffixe_array.append(new_tup)
            idx+=1
        suffixe_array.sort()
        partial_suffixe_array = []
        if partial == True:
            idx = 0
            for elem in suffixe_array:
                if elem[1]%k == 0:
                    # We got a multiple of k, inserting this elem in the array.
                    # idx represent position from the original suffixe array, and elem[1] is the position
                    # of the suffixe in the orignal text.
                    partial_suffixe_array.append((idx,elem[1]))
                idx+=1
            return partial_suffixe_array
        return suffixe_array

    def sort_list_suffixe_array(self, new_list):
        """
            The function will sort all the tuple present in the list.
        """
        return 1

    def print_idx_suffixe_array(self, list, partial=False):
        """
            The function will only print the location of
        """
        for elem in list:
            if partial == False:
                print(str(elem[1])+',',)
            else:
                print(str(elem[0])+','+str(elem[1]))

    def last_to_first(self, matrix_bw):
        """
            Function used for BW_Matching. Will the return the index on the first column,
                A1 ... B1
                B1 ... A1
                B2 ... B2
                A2 ... B1
            -> By instance , last_to_first(0) = 1 , because char on line 0 of last column is B1, and this B1
                can be find in the first column at index 1
        """
        list_last_column = matrix_bw[len(matrix_bw[0])-1]
        list_first_column = matrix_bw[0]
        list_last_to_first = []
        count_A = 1
        count_T = 1
        count_C = 1
        count_G = 1
        for char in list_last_column:
            if char == 'A':
                idx_mapping_to_first_column = list_first_column.index(( char, count_A))
                count_A+=1
            elif char == 'T':
                idx_mapping_to_first_column = list_first_column.index(( char, count_T))
                count_T+=1
            elif char == 'C':
                idx_mapping_to_first_column = list_first_column.index(( char, count_C))
                count_C+=1
            elif char == 'G':
                idx_mapping_to_first_column = list_first_column.index(( char, count_G))
                count_G+=1
            list_last_to_first.append(idx_mapping_to_first_column)
        return list_last_to_first

    def transform_list_to_string(self, list):
        """
            Self-explanatory.
        """
        text = ''
        for x in list:
            text = text + str(x)
        return text

    def bw_matching(self, first_column, last_column, pattern, last_to_first, better_bw_matching=False):
        """
            Implement the algo defined on coursera.org
            Number of pattern present in a big string.
            last_column should be first converted into a string, because findint the first and last occurence
            is easier with a string.
        """
        if better_bw_matching == True:
            self.build_matrix_BW_count(last_column)
            first_occurence = self.first_occurence_first_column( first_column)
        last_column_text = self.transform_list_to_string(last_column)
        top = 0
        # To get the relative position [ comparing with orig text ] of first occurence of symbol.
        total_count_top = 0
        bottom = len(last_column_text)
        while top <= bottom:
            if len(pattern) > 0:
                idx_last_letter = len(pattern)-1
                symbol = pattern[idx_last_letter]
                #pattern = pattern[0:idx_last_letter]
                pattern = buffer(pattern, 0, idx_last_letter)
                if symbol in last_column_text[top:bottom]:
                    if better_bw_matching == False:
                        # The position should always be in comparison to the all list, that's why we add
                        # top , at start it is 0 -> no prob !
                        #top_index = ( last_column_text.find(symbol, top, bottom) ) + total_count_top
                        top_index = ( last_column_text.find(symbol, top, bottom) )
                        bottom_index = last_column_text.rfind(symbol, top, bottom)
                        total_count_top = total_count_top + top_index
                        top = last_to_first[top_index]
                        # we put a +1, otherwise, the search does not take into account this last character for the
                        # next search which will be [top , bottom ] , bottom previously represented the last
                        # character, but if a search with slicing, it does not take it !
                        bottom = last_to_first[bottom_index] + 1
                    else:
                        # implementing with the method of the better BW matching.
                        idx_symbol = self.get_idx_character( symbol)
                        top = first_occurence[idx_symbol] + 1
                        bottom = 1
                else:
                    return 0
            else:
                # -1 removed, does not give the good result with that, except with 0
                #return bottom-top+1
                return bottom-top

    def get_idx_character(self, character):
        """
            Function just used to translate a character to a numeric value,
            which correspond to an index in a matrices.
        """
        if character == '$':
            return 0
        elif character == 'A':
            return 1
        elif character == 'T':
            return 2
        elif character == 'C':
            return 3
        elif character == 'G':
            return 4

    def increment_count_char( self, character):
        """
            Function use to update the count of the character seen in a list.
            We define idx column for character as : $=0, A=1, T=2, C=3, G=4
        """
        idx = self.get_idx_character(character)
        current_count_char = self.count_char[idx]
        current_count_char+=1
        self.count_char[idx] = current_count_char
        return self.count_char[idx]

    def get_count_character(self, idx):
        """
            Function used to get the current count of character seen so far in a list.
            possible character are : $,A,T,C,G -> 5
        """
        return self.count_char[idx]

    def build_matrix_BW_count(self, list_last_column_matrix):
        """
            Function used in BW matching to be faster; as said in coursera.org.
            Instead of searchning each time for a symbol in the last column, we
            construct a matrix. count will return the number of symbol so far seen until
             idx in the last column.
        """
        # 5 for the matrix, because : $ A T C G
        # len + 1 , because at position 0 nothing is seen -> all is zero.
        self.bw_count = numpy.zeros( (len(list_last_column_matrix)+1,5) )
        # Beginning filling the matrix from idx 1 because on first row, all should be 0.
        idx = 1
        # we define idx column for character as : $=0, A=1, T=2, C=3, G=4
        for char in list_last_column_matrix:
            idx_column = self.get_idx_character(char)
            new_count_char = self.increment_count_char( char)
            self.bw_count[idx][idx_column] = new_count_char
            idx+=1

    def first_occurence_first_column(self, list):
        """
            Create list with the position of first occurence of a char in a list.
            Possible char are : $ A T C G
            # 5 for the list , because : $ A T C G , respectively idx are 0 1 2 3 4
        """
        list_first_occurence = []
        list_nucleotide.append('$')
        list_nucleotide.extend(['A','T','C','G'])
        for nucleotide in list_nucleotide:
            # We should pay attention to the fact, that may be one nucleotide is not present ( may be... )
            # ATTENTION !
            idx_char = self.get_idx_character(nucleotide)
            list_first_occurence[idx_char] = list.index(nucleotide)
        return list_first_occurence

    def build_checkpoint_array(self, last_column, k):
        """
            We will be a matrix containing only the idx of the multiple K.
            Each row will contains the number of character so far seen in the list last_column.
            Each column is a count of these char, idx respectively : $,A,T,C,G -> 0,1,2,3,4
        """
        checkpoint_array = numpy.zeros( (len(last_column),5) )
        # The first row of the matrix should be filled with 0 as said in coursera.org [ stepic.org ]
        # That's why we start with 1
        idx = 1
        # Re-initializing the list with the total count of each char so far seen in the list.
        # Because this list is also used for another function.
        self.count_char = [0,0,0,0,0]
        for char in last_column:
            idx_column = self.get_idx_character(char)
            current_count = self.increment_count_char(char)
            checkpoint_array[idx][idx_column] = current_count
            idx+=1