# Class to create a codon table for humans

# data from http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606&aa=1&style=N retrieved on 11/20/2014

# accessor methods to return avg codon frequency for a sequence or to find the stretch of N codons with minimum frequency, and what that frequency is


class HumanCodonTable:

    def __init__(self):
        # define the codon table here.

        self.codonTable = {'TTT': .46,'TCT':0.19,'TCC':0.22,'TCA':0.15,'TCG':0.05,'CCT':0.29,'CCC':0.32,'CCA':0.28,'CCG':0.11,'ACT':0.25,'ACC':0.36,'ACA':0.28,'ACG':0.11\
                               ,'GCT':0.27,'GCC':0.4,'GCA':0.23,'GCG':0.11,'TAT':0.44,'TAC':0.56,'TAA':0.3,'TAG':0.24,'CAT':0.42,'CAC':0.58,'CAA':0.27,'CAG':0.73,'AAT':0.47\
                               ,'AAC':0.53,'AAA':0.43,'AAG':0.57,'GAT':0.46,'GAC':0.54,'GAA':0.42,'GAG':0.58,'TGT':0.46,'TGC':0.54,'TGA':0.47,'TGG':1,'CGT':0.08,'CGC':0.18\
                               ,'CGA':0.11,'CGG':0.2,'AGT':0.15,'AGC':0.24,'AGA':0.21,'AGG':0.21,'GGT':0.16,'GGC':0.34,'GGA':0.25,'GGG':0.25,'TTT':0.46,'TTC':0.54,'TTA':0.08\
                               ,'TTG':0.13,'CTT':0.13,'CTC':0.2,'CTA':0.07,'CTG':0.4,'ATT':0.36,'ATC':0.47,'ATA':0.17,'ATG':1,'GTT':0.18,'GTC':0.24,'GTA':0.12,'GTG':0.46}
    # calculate the average codon usage frequency across the entire input sequence 

    def averageCodonFreq(self, sequence): 
        numCodons = len(sequence)/3
        return sum([self.codonTable[sequence[i*3:i*3+3]] for i in range(numCodons)])/numCodons


    # note - codonWindowSize here is interpreted as a number of codons, hence multiplied by 3 inside. 
    def minCodonFreq(self, sequence, codonWindowSize):
        numCodons = len(sequence)/3
        if (numCodons < codonWindowSize): 
            codonWindowSize = numCodons
        allWindows = [sequence[i*3:(i+codonWindowSize)*3] for i in range(numCodons - codonWindowSize + 1)]
        allFrequencies = map(self.averageCodonFreq, allWindows)
        minFrequency = min(allFrequencies) 

        return minFrequency 
