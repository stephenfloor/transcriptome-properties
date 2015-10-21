# converts between two annotation sets.  First column of input is assumed to be the "from" annotations, and second column is the "to" annotations

# i.e. 

# ENST00000517143 NR_046932

# will convert from ensembl to refseq 


from collections import defaultdict 

class AnnotationConverter:
    def __init__(self, annotConvFile):
        
        self.conversion = defaultdict(list)
        
        for line in annotConvFile:
            line = line.split()
            self.conversion[line[0]].append(line[1])
            
    def convert(self, id):
        #print "AnnotConv: %s maps to %s" % (id, self.conversion[id])
        return self.conversion[id]

