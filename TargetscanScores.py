# TargetscanScores.py 

# class to read a targetscan Summary_Counts.txt type file and parse it into gene-specific information with reasonable accessor methods 

from collections import defaultdict

class TargetscanScores: 
    def __init__(self, targetscanScoreFile): 
        
        # compute: 
        # 1) sum of context+ scores 
        # 2) min context+ score 
        # 3) # of sites
        # 4) sum of conserved site context+ scores
        # 5) min conserved context+ score
        # 6) # of conserved sites 


        self.refSeqToScores = defaultdict(list)  # stores a list of all miRNA sites per refseq ID 
        self.refSeqToScoreSum = defaultdict(float) # stores the sum of context+ scores per refseq ID 
        self.refSeqToMinScore = defaultdict(float) # stores the min context+ score per refseq ID 
        self.refSeqToNumSites = defaultdict(int) # stores the number of sites per refseq ID 
        self.refSeqToCnsvScores = defaultdict(list) # stores a list of all conserved miRNA sites per refseq ID 
        self.refSeqToCnsvScoreSum = defaultdict(float) # stores a list of the sum of scores for conserved site per refseq ID 
        self.refSeqToCnsvMinScore = defaultdict(float) # stores the min conserved score per refseq ID 
        self.refSeqToCnsvNumSites = defaultdict(int) # stores the number of conserved sites
        self.refSeqToNoncnsvNumSites = defaultdict(int) # stores the number of nonconserved sites
        self.refSeqToNoncnsvScores = defaultdict(list) # stores a list of all conserved miRNA sites per refseq ID 
        self.refSeqToNoncnsvScoreSum = defaultdict(float) # stores a list of the sum of scores for conserved site per refseq ID 
        self.refSeqToNoncnsvMinScore = defaultdict(float) # stores the min conserved score per refseq ID 
        self.refSeqToCnsv8merSites = defaultdict(int) # stores the number of conserved 8mer sites
        self.refSeqToCnsv7merm8Sites = defaultdict(int) # cnsv 7mer-m8 sites
        self.refSeqToCnsv7mer1aSites = defaultdict(int) # cnsv 7mer-1a sites 
        self.refSeqToNoncnsv8merSites = defaultdict(int) # stores the number of nonconserved 8mer sites
        self.refSeqToNoncnsv7merm8Sites = defaultdict(int) # noncnsv 7mer-m8 sites
        self.refSeqToNoncnsv7mer1aSites = defaultdict(int) # noncnsv 7mer-1a sites 

        
        
#0 Transcript ID   1 Gene Symbol     2 miRNA family    3 Species ID      4 Total num conserved sites       5 Number of conserved 8mer sites  6 Number of conserved 7mer-m8 sites       7 Number of conserved 7mer-1a sites
#       8 Total num nonconserved sites    9 Number of nonconserved 8mer sites      10 Number of nonconserved 7mer-m8 sites    11 Number of nonconserved 7mer-1a sites    12 Representative miRNA    13 Total context score     14 Aggregate PCT
#NM_000014       A2M     AAAAGUG 9606    0       0       0       0       1       0       1       0       hsa-miR-548t    -0.128  NULL

        # skip the header
        targetscanScoreFile.readline()


        for line in targetscanScoreFile:
            if (not line.strip()):
                continue
            
            line = line.split() 
            if (line[13] == "NULL"):
                continue

            isConserved = int(line[4])
            
            self.refSeqToScores[line[0]].append( (line[12], float(line[13]) ) )
            self.refSeqToScoreSum[line[0]] += float(line[13])
            self.refSeqToNumSites[line[0]] += int(line[4]) + int(line[8])
            
            if (isConserved): 
                self.refSeqToCnsvScores[line[0]].append( (line[12], float(line[13]) ) )
                self.refSeqToCnsvScoreSum[line[0]] += float(line[13])
                self.refSeqToCnsv8merSites[line[0]] += int(line[5])
                self.refSeqToCnsv7merm8Sites[line[0]] += int(line[6])
                self.refSeqToCnsv7mer1aSites[line[0]] += int(line[7]) # cnsv 7mer-1a sites 
                self.refSeqToCnsvNumSites[line[0]] += int(line[4])

            else: 
                self.refSeqToNoncnsvScores[line[0]].append( (line[12], float(line[13]) ) )
                self.refSeqToNoncnsvScoreSum[line[0]] += float(line[13])
                self.refSeqToNoncnsv8merSites[line[0]] += int(line[9])
                self.refSeqToNoncnsv7merm8Sites[line[0]] += int(line[10])
                self.refSeqToNoncnsv7mer1aSites[line[0]] += int(line[11]) # cnsv 7mer-1a sites 
                self.refSeqToNoncnsvNumSites[line[0]] += int(line[8])
                
                
#0 Transcript ID   1 Gene Symbol     2 miRNA family    3 Species ID      4 Total num conserved sites       5 Number of conserved 8mer sites  6 Number of conserved 7mer-m8 sites       7 Number of conserved 7mer-1a sites
#       8 Total num nonconserved sites    9 Number of nonconserved 8mer sites      10 Number of nonconserved 7mer-m8 sites    11 Number of nonconserved 7mer-1a sites    12 Representative miRNA    13 Total context score     14 Aggregate PCT

        # compute meta properties here (min score, etc) 
        for key,val in self.refSeqToScores.iteritems():
            self.refSeqToMinScore[key] = min([score[1] for score in val])

        for key,val in self.refSeqToCnsvScores.iteritems():
            self.refSeqToCnsvMinScore[key] = min([score[1] for score in val])

        for key,val in self.refSeqToNoncnsvScores.iteritems():
            self.refSeqToNoncnsvMinScore[key] = min([score[1] for score in val])


    def getScores(self, refSeqGene):
        return self.refSeqToScores[refSeqGene]
    
    def getScoreSum(self, refSeqGene): 
        return self.refSeqToScoreSum[refSeqGene]
    
    def getMinScore(self, refSeqGene): 
        return self.refSeqToMinScore[refSeqGene]

    def getNumSites(self, refSeqGene):
        return self.refSeqToNumSites[refSeqGene]
    
    def getCnsvScores(self, refSeqGene):
        return self.refSeqToCnsvScores[refSeqGene]
    
    def getCnsvScoreSum(self, refSeqGene):
        return self.refSeqToCnsvScoreSum[refSeqGene]

    def getCnsvMinScore(self, refSeqGene):
        return self.refSeqToCnsvMinScore[refSeqGene]
    
    def getCnsvNumSites(self, refSeqGene):
        return self.refSeqToCnsvNumSites[refSeqGene]
    
    def getNoncnsvNumSites(self, refSeqGene):
        return self.refSeqToNoncnsvNumSites[refSeqGene]

    def getNoncnsvScores(self, refSeqGene):
        return self.refSeqToNoncnsvScores[refSeqGene]

    def getNoncnsvScoreSum(self, refSeqGene):
        return self.refSeqToNoncnsvScoreSum[refSeqGene]

    def getNoncnsvMinScore(self, refSeqGene):
        return self.refSeqToNoncnsvMinScore[refSeqGene]

    def getCnsv8merSites(self, refSeqGene):
        return self.refSeqToCnsv8merSites[refSeqGene]

    def getCnsv7merm8Sites(self, refSeqGene):
        return self.refSeqToCnsv7merm8Sites[refSeqGene]

    def getCnsv7mer1aSites(self, refSeqGene):
        return self.refSeqToCnsv7mer1aSites[refSeqGene]

    def getNoncnsv8merSites(self, refSeqGene):
        return self.refSeqToCnsv8merSites[refSeqGene]

    def getNoncnsv7merm8Sites(self, refSeqGene):
        return self.refSeqToCnsv7merm8Sites[refSeqGene]

    def getNoncnsv7mer1aSites(self, refSeqGene):
        return self.refSeqToCnsv7mer1aSites[refSeqGene]

            
#         print self.refSeqToScores 
#         print self.refSeqToScoreSum
#         print self.refSeqToMinScore
#         print "self.refSeqToNumSites"
#         print self.refSeqToNumSites
#         print "self.refSeqToCnsvScores"
#         print self.refSeqToCnsvScores
#         print "self.refSeqToCnsvScoreSum"
#         print self.refSeqToCnsvScoreSum
#         print "self.refSeqToCnsvMinScore"
#         print self.refSeqToCnsvMinScore
#         print "self.refSeqToNumCnsvSites"
#         print self.refSeqToCnsvNumSites
#         print "self.refSeqToNoncnsvNumSites"
#         print self.refSeqToNoncnsvNumSites
#         print "self.refSeqToNoncnsvScores"
#         print self.refSeqToNoncnsvScores
#         print "self.refSeqToNoncnsvScoreSum"
#         print self.refSeqToNoncnsvScoreSum
#         print "self.refSeqToNoncnsvMinScore"
#         print self.refSeqToNoncnsvMinScore
#         print "self.refSeqToCnsv8merSites"
#         print self.refSeqToCnsv8merSites
#         print "self.refSeqToCnsv7merm8Sites"
#         print self.refSeqToCnsv7merm8Sites
#         print "self.refSeqToCnsv7mer1aSites"
#         print self.refSeqToCnsv7mer1aSites
#         print "self.refSeqToNoncnsv8merSites"
#         print self.refSeqToNoncnsv8merSites
#         print "self.refSeqToNoncnsv7merm8Sites"
#         print self.refSeqToNoncnsv7merm8Sites
#         print "self.refSeqToNoncnsv7mer1aSites"
#         print self.refSeqToNoncnsv7mer1aSites

            
