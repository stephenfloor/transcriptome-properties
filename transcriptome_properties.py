#!/usr/bin/env python


# transcriptome_properties.py

# Calculate (ir)relevant properties of an input transcriptome (or subset thereof) specified in BED format

# TODO 

# make sure # of header coluns = # of output columns for all files 
# search for todo throughout 

# Properties to calculate:

# GC content (done)
# length (done)
# exonct (done) 

# optional (region-specific):

# - 5' UTR specific
#  - structure at 5' end (for 5' UTR, specifically affects 43S loading)
#  - min deltaG over sliding window of N nts (default 75)
#  - 5' TOP or etc?  pyrimidine content of 5' proximal?  sabatini defined this as a 5nt pyrimidine stretch within 15nt of TSS. (not implemented) 

# start codon specific 
#  - kozak context (done) 
#  - number of uORFs (done)
#  - uorf overlap with start codon (done)
#  - code all four possibilities for each ATG-proximal base (not implemented) 
#  - start codon itself - is it always ATG?  (done)

# - CDS specific
#  - % rare codons 
#  - signal peptide?  (not implemented; just download list from http://www.signalpeptide.de/index.php  ? ) 
#  - min deltaG over sliding window of N nts


# - 3' UTR specific
#  - miRNA binding sites (targetscan)
#  - AU-rich elements
#  - min deltaG over sliding window of N nts 

import sys,os,argparse,subprocess,re,multiprocessing,threading,time
import HumanCodonTable
import AnnotationConverter
import TargetscanScores
from SNFUtils import *
from collections import defaultdict
from threading import Lock     

# gets the folding energy of a sequence. returns a tuple containing energies (MFE, centroid, MEA) 
# RNAfold courtesy viennaRNA 
def get_energy(sequence): 
    if len(sequence) < 1: 
        return (0, 0, 0)

    lines = stdout_from_command("echo %s | RNAfold --noPS --MEA -p" % sequence)

    # skip the first line of lines - just the seq 
    next(lines)

    #second line has the MFE 
    MFE = float(re.sub('[()]','', next(lines).decode('utf8').split()[-1]))
        
    #third line has the ensemble energy - discarding for now
    next(lines)
    
    #fourth line has the centroid energy 
    centroid = float(re.sub('[{}]','', next(lines).decode('utf8').split()[-2]))
    
    #fifth line has the MEA 
    MEA = float(re.sub('[{}]','', next(lines).decode('utf8').split()[-2]))

    return (MFE, centroid, MEA)

# gets the min folding energy for a sliding window of windowsize across the sequence
# uses multiprocessing via the procpool to accelerate. - removed in favor of global threading
#def get_sliding_energy(sequence, windowsize, procpool):
def get_sliding_energy(sequence, windowsize):
    if len(sequence) < 1: 
        return (0, 0, 0)
    all75s = [sequence[i:i+windowsize] for i in range(len(sequence) - windowsize)]
    #nrg = procpool.map(get_energy, all75s)
    nrg = map(get_energy, all75s)
    minMFE = min([tup[0] for tup in nrg])
    mincentroid = min([tup[1] for tup in nrg])
    minMEA = min([tup[2] for tup in nrg])
    
    return (minMFE, mincentroid, minMEA)


# gets the folding energy (MFE) of a sequence using RNALfold
def get_rnalfold_energy(sequence, windowsize): 
    if len(sequence) < 1: 
        return 0

    lines = stdout_from_command("echo %s | RNALfold -L %d" % (sequence, windowsize))

    # example output: .(((.(((..((((.....((((....((.......))....)))).)))).............))).))). ( -9.00)   44
    minMFE = min([float(line.split(' (')[1][:6]) if line[0] in ".(" else "" for line in lines])

    if not minMFE:
        minMFE = 0

    return float(minMFE)
                                

parser = argparse.ArgumentParser(description="Calculate transcriptome-wide properties")

globalargs = parser.add_argument_group(title="Global arguments")
globalargs.add_argument("-i", "--input", help="The transcriptome region (BED format)", required=True)
globalargs.add_argument("-g", "--genome", help="The genome for the input transcriptome", required=True)
globalargs.add_argument("--gc", help="Calculate GC content", action="store_true")
globalargs.add_argument("--length", help="Calculate length", action="store_true")
globalargs.add_argument("--exonct", help="Count # of exons", action="store_true") 
globalargs.add_argument("--nt", help="Number of threads (default is 8 or 4 for lfold)", default=8, type=int)
globalargs.add_argument("-o", "--output", help="Output basename (e.g. CDS)", default="tx_props")
globalargs.add_argument("--window", help="Window size for sliding window calculations (default 75)", default=75, type=int)
globalargs.add_argument("--convtorefseq", help="Filename to convert input annotations to refseq (for targetscan; e.g. knownToRefSeq.txt)", default=False)
globalargs.add_argument("--targetscanfile", help="Filename of targetscan scores (e.g. Summary_Counts.txt)", default=False)
globalargs.add_argument("--deltag", help="Calculate min deltaG in sliding window of size --window over region", action="store_true")
globalargs.add_argument("--lfold", help="Use RNALfold to calculate MFE rather than RNAfold (faster but does not compute centroid,MEA)", action="store_true")

utr5args = parser.add_argument_group(title="5' UTR specific arguments")
utr5args.add_argument("--cap-structure", help="Calculate structure at the 5' end", action="store_true")
#utr5args.add_argument("--deltag-5utr", help="Calculate min deltaG in 5' UTR in sliding window of size --window", action="store_true")

startargs = parser.add_argument_group(title="Start-codon-specific arguments")
startargs.add_argument("--kozak", help="Calculate Kozak context score", action="store_true")
startargs.add_argument("--uorf-count", help="Calculate number of 5' UTR uORFs (starting with [ACT]TG)", action="store_true")
startargs.add_argument("--uorf-overlap", help="Overlap of uORF with start codon (implies --uorf-count)", action="store_true")
startargs.add_argument("--start-codon", help="Record the start codon used (ATG or other)", action="store_true") 

cdsargs = parser.add_argument_group(title="CDS-specific arguments") 
cdsargs.add_argument("--rare-codons", help="Calculate codon usage properties", action="store_true") 
#cdsargs.add_argument("--deltag-cds", help="Calculate min deltaG in CDS in sliding window of size --window", action="store_true") 

utr3args = parser.add_argument_group(title="3' UTR specific arguments") 
utr3args.add_argument("--mirna-sites", help="Compile miRNA binding site info from targetscan", action="store_true")
utr3args.add_argument("--au-elements", help="Count number of AU-rich elements in the 3' UTR", action="store_true")
#utr3args.add_argument("--deltag-3utr", help="Calculate min deltaG in 3' UTR in sliding window of size --window", action="store_true") 

args = parser.parse_args()

print("-----------------------------------")
print("|   transcriptome_properties.py   |")
print("|    compute props of a txome     |")
print("|      run with -h for help       |")
print("|          snf 11/2014            |")
print("-----------------------------------\n\n")


print ("Arguments:")

print (args)

# get a streaming input of all the regions from bedtools

print( "\n\nReading from regions using bedtools getfasta: ")

# bedtools getfasta flags: 
#  -name  preserve region name
#  -tab   tab-delimited output 
#  -s     reverse complement - strand regions 
#  -split concat all blocks in the input bed (i.e. splice together exons) 

cmd = "bedtools getfasta -fi %s -bed %s -name -tab -s -fo - -split 2> %s_stderr.log" % (args.genome, args.input, args.output)
print( cmd + "\n")

# sanity check inputs here

prompted = False

if (args.lfold):
    print ("Using RNALfold to calculate just MFE\n")
    if (args.nt != 4):
        print ("Warning: highest performance during benchmarks with --lfold achieved with 4 threads\n")
    

if (args.cap_structure):
    prompt("Using 5' UTR specific options, so input should be 5' UTR regions!")
    prompted = True

if (args.kozak or args.uorf_overlap or args.uorf_count or args.start_codon):
    if (prompted):
        sys.exit("FATAL: multiple region-specific options for disparate regions selected.")
    prompt("Using start codon specific options, so input should be start codon regions (5utr_start)!") 
    prompted = True 

if (args.rare_codons):
    if (prompted):
        sys.exit("FATAL: multiple region-specific options for disparate regions selected.")
    prompt("Using CDS specific options, so input should be CDS regions!") 
    prompted = True 

if (args.mirna_sites):
    if (prompted):
        sys.exit("FATAL: multiple region-specific options for disparate regions selected.")
    prompt("Using 3' UTR specific options, so input should be 3' UTR regions!") 
    prompted = True 

if (args.mirna_sites): 
    if (not args.convtorefseq):
        sys.exit("FATAL: argument --convtorefseq <file> required for --mirna-sites")
    if (not args.targetscanfile):
        sys.exit("FATAL: argument --targetscanfile <file> required for --mirna-sites") 

# open relevant output files 

if (args.gc):
    gc_outfile = safe_open_file("%s_gc_content.csv" % args.output)
    gc_outfile.write("transcriptID,gc_content\n")

if (args.length):
    len_outfile = safe_open_file("%s_length.csv" % args.output)
    len_outfile.write("transcriptID,length\n")

txid_to_exonct = defaultdict(int)

if (args.exonct):
    exonct_outfile = safe_open_file("%s_exonct.csv" % args.output)
    exonct_outfile.write("transcriptID,exonct\n")
    print( "Building exon count dictionary...")
    # count the # of exons real quick and throw it into a dictionary
    # col 9 is exonct, col 3 is txid 
    with open(args.input, "rt") as inp:
        for line in inp:
            line = line.split()
            txid_to_exonct[line[3]] = line[9]
    print ("Built.\n")

if (args.cap_structure):
    cap_structure_outfile = safe_open_file("%s_cap_structure.csv" % args.output) 
    cap_structure_outfile.write("transcriptID,MFE,centroid,MEA\n")

if (args.deltag):
    if (args.lfold): 
        deltag_outfile = safe_open_file("%s_deltag_lfold.csv" % args.output)
        deltag_outfile.write("transcriptID,MFE_min\n")
    else:
        deltag_outfile = safe_open_file("%s_deltag.csv" % args.output)
        deltag_outfile.write("transcriptID,MFE_min,centroid_min,MEA_min\n")

if (args.uorf_count or args.uorf_overlap): # output the count even if just overlap is supplied; have to calculate this anyway
    uorf_count_outfile = safe_open_file("%s_uorf_count.csv" % args.output) 
    uorf_count_outfile.write("transcriptID,uORF_count\n") 

if (args.uorf_overlap):
    uorf_overlap_outfile = safe_open_file("%s_uorf_overlap.csv" % args.output)
    uorf_overlap_outfile.write("transcriptID,uORF_overlap\n")

if (args.kozak):
    kozak_outfile = safe_open_file("%s_kozak.csv" % args.output) 
    kozak_outfile.write("transcriptID,kozak_score\n")

if (args.start_codon): 
    start_codon_outfile = safe_open_file("%s_start_codon.csv" % args.output)
    start_codon_outfile.write("transcriptID,start_codon,start_codon_numeric\n")

if (args.rare_codons):
    rare_codons_outfile = safe_open_file("%s_rare_codons.csv" % args.output)
    rare_codons_outfile.write("transcriptID,avg_codon_freq,min_codon_freq\n") 

if (args.mirna_sites):
    mirna_sites_outfile = safe_open_file("%s_mirna_count.csv" % args.output)
    mirna_sites_outfile.write("transcriptID,num_miRNA_sites,miRNA_score_sum,miRNA_score_min,num_cnsv_sites,cnsv_score_sum,cnsv_score_min,num_noncnsv_sites,noncnsv_score_sum,noncnsv_score_min\n")
    mirna_density_outfile = safe_open_file("%s_mirna_density.csv" % args.output)
    mirna_density_outfile.write("transcriptID,miRNA_density,miRNA_score_density,cnsv_site_density,cnsv_score_density,noncnsv_site_density,noncnsv_score_density\n")

    try: 
        convFile = open(args.convtorefseq, 'r')
        conv2rs = AnnotationConverter.AnnotationConverter(convFile)
    except: 
        sys.exit("Cannot open file %s for annotation conversion to refseq" % args.convtorefseq)
    try:
        targetscanScoreFile = open(args.targetscanfile, 'r')
    except: 
        sys.exit("Cannot open file %s for reading targetscan scores" % args.targetscanfile)
    tsscores = TargetscanScores.TargetscanScores(targetscanScoreFile)


if (args.au_elements): 
    au_elements_outfile = safe_open_file("%s_au_elements.csv" % args.output)
    au_elements_outfile.write("transcriptID,pentamer_count,au_element_count,au_element_frac,max_au_length\n")

#pool = multiprocessing.Pool(processes=args.nt)


def process_line(line, lock):

    line = line.split()

    transcriptID = line[0]
    #print "Processing tx %s" % transcriptID

    sequence = line[1].upper()

    seqlen = len(sequence)

    if (args.gc):

        numN = sequence.count("N")

        if ("N" in sequence): 
            #print "Warning: stripping N nucleotides from transcript %s" % (transcriptID) 
            sequence = sequence.strip("N") 
            #seqlen = len(sequence) 

        numG = sequence.count("G")
        numC = sequence.count("C")
        numA = sequence.count("A")
        numT = sequence.count("T")

        #if (numG + numC + numA + numT != seqlen):  # there must be an internal N; try to continue somewhat gracefully by adjusting the sequence length
            #seqlen -= numN
            #print "Warning: %d internal N nucleotides ignored when calculating GC content for transcript %s" % (numN, transcriptID) 
            #sys.exit("FATAL: sum of ACTG (%d) != length of sequence (%d) for sequence %s" % (numG+numC+numA+numT, seqlen, sequence))
            
        percent_gc = (numG + numC) / float(seqlen - numN) if seqlen - numN > 0 else 0

        '''acquire a thread lock'''
        lock.acquire()
        gc_outfile.write("%s,%3.2f\n" % (transcriptID, percent_gc))

        '''release the lock'''
        lock.release()


    if (args.length):
        lock.acquire()
        len_outfile.write("%s,%d\n" % (transcriptID, seqlen))
        lock.release()

    if (args.exonct):
        lock.acquire()
        exonct_outfile.write("%s,%s\n" % (transcriptID, txid_to_exonct[transcriptID.split("::")[0]]))
        lock.release()
        
    if (args.cap_structure): 
        # calculate the deltaG of the 50nt after the 5' end here, or if the 5' UTR is less than 50nt just calculate the deltaG of the whole thing.
        # using viennaRNA 
        if (seqlen < 50):
            nrg = get_energy(sequence)
        else:
            nrg = get_energy(sequence[0:50])

        lock.acquire()
        cap_structure_outfile.write("%s,%.1f,%.1f,%.1f\n" % (transcriptID, nrg[0],nrg[1],nrg[2]))
        lock.release()
        

    if (args.deltag):
        # calculate the min deltag of a sliding window of size args.window across the region 
        if (args.lfold):
            lock.acquire()
            deltag_outfile.write("%s,%.1f\n" % (transcriptID, get_rnalfold_energy(sequence, args.window)))
            lock.release()
        else: 
            if (seqlen <= args.window):
                nrg = get_energy(sequence)
                lock.acquire()
                deltag_outfile.write("%s,%.1f,%.1f,%.1f\n" % (transcriptID, nrg[0],nrg[1],nrg[2]))
                lock.release()
            else: 
                #nrg = get_sliding_energy(sequence, args.window, pool)
                nrg = get_sliding_energy(sequence, args.window)
                lock.acquire()
                deltag_outfile.write("%s,%.1f,%.1f,%.1f\n" % (transcriptID, nrg[0],nrg[1],nrg[2]))
                lock.release()
        
    if (args.kozak): 
        # score the Kozak context.  This is G c c A/G c c atg G.  The most important nts are +4, -3 and -6.  Scoring these as +3 and the others as +1. Max score = 13
        # need to implement this.  this is actually a "start codon" specific parm because it involves overlap between the 5' UTR and CDS. likewise for uorf-overlap 
        
        if (len(sequence) > 32): 
            # start codon is at -28 through -26 in the 5utr_start files, so kozak starts at -34 and ends at -24 ignoring -28 through -26  
            kozakScore = sum ((sequence[-32] == "C", sequence[-31] == "C", sequence[-29] == "C", sequence[-28] == "C"))
            kozakScore += 3 * sum((sequence[-33] == "G", sequence[-30] == "A" or sequence[-30] == "G", sequence[-24] == "G"))
            
        else: 
            # something with a UTR that is less than 6nt is not likely to follow kozak behavior, so just setting the score to -1 to indicate oddball status 
            kozakScore = -1
    
        lock.acquire()
        kozak_outfile.write("%s,%d\n" % (transcriptID, kozakScore))
        lock.release()

    if (args.uorf_count or args.uorf_overlap): 
        # count the number of uORFs (min length 10 codons) in the sequence. 
        # !!! IT IS ASSUMED that the sequence is the 5' UTR ending with ATG([N]{27}) as in the 5utr_start bed files, otherwise this makes no sense!!!!!!!!!!
        
        # uORF is defined as orf >= 10 codons that starts with ATG/CTG/GTG as per Ingolia 2011 and either ends with a stop codon or the end of the UTR which implies it continues into the proper CDS 
        # I am not sure how good these criteria really are, and they predict massive numbers of uORFs. an alternative approach would be to actually measure the translated uORFs with harringtonine in these cells... 
        

        # re to find all uORFs including near-cognate matches CTG/TTG 
        orf_pattern = re.compile(r'(?=([ACG]TG(?:...)*?)(?=TAG|TGA|TAA|.{0,2}$))')  # EOL counts because that means the uORF has read at least 27nt into the coding sequence
        uorfs = [orf for orf in orf_pattern.findall(sequence) if len(orf) >= 30]

        # trim uorfs - method one; use string comparisons 

        uorfs_trimmed=[]
        for orf in uorfs:
            found = False 
            if (orf not in uorfs_trimmed):
                for orf_cmp in uorfs_trimmed:
                    if (orf in orf_cmp):
                        found = True 
                if (not found):
                    uorfs_trimmed.append(orf)
            
        uorfs = uorfs_trimmed

        # nice up the output
        lock.acquire()
        uorf_count_outfile.write("%s,%d\n" % (transcriptID, len(uorfs)))
        lock.release()

        if (args.uorf_overlap): 
            
            starts = [sequence.find(orf) for orf in uorfs]
            ends = [starts[i] + len(uorfs[i]) for i in range(len(uorfs))]
            overlaps = sum([seqlen - ends[i] < 28 for i in range(len(ends))]) # 28 here comes from the construction of 5utr_start bed files, which have 27nt after the start.

            #print sequence 
            #print starts, ends, overlaps
            
            lock.acquire()
            uorf_overlap_outfile.write("%s,%d\n" % (transcriptID, overlaps))
            lock.release()


    if (args.start_codon): 
        # record what the start codon is, along with a numeric representation (0 is ATG, 1 is 1 mismatch, 2 is 2 mismatches, 3 is 3 mismatches)
        # INPUT MUST BE 5utr_start for these calculations, in which the start codon is from 
        

        startCodon = sequence[-27:-24]

        mismatches = sum((startCodon[0] != "A", startCodon[1] != "T", startCodon[2] != "G"))

        #print sequence, startCodon, mismatches

        lock.acquire()
        start_codon_outfile.write("%s,%s,%d\n" % (transcriptID, startCodon, mismatches))
        lock.release()


    if (args.rare_codons):
        codonTable = HumanCodonTable.HumanCodonTable()
        
        lock.acquire()
        rare_codons_outfile.write("%s,%3.2f,%3.2f\n" % (transcriptID, codonTable.averageCodonFreq(sequence), codonTable.minCodonFreq(sequence, 5)))  # 5 here is the window size in codons, so 15 nucleotides 
        lock.release()



    if (args.mirna_sites):

        refseqIDs = conv2rs.convert(transcriptID[:-5])

        if (refseqIDs): 

            nSites = 0 
            scoreSum = 0
            minScore = 999999999
            nCnsvSites = 0 
            cnsvScoreSum = 0
            cnsvMinScore = 9999999999
            nNoncnsvSites = 0 
            noncnsvScoreSum = 0
            noncnsvMinScore = 9999999999

            
            for id in refseqIDs: 
                nSites += tsscores.getNumSites(id)
                scoreSum += tsscores.getScoreSum(id) 
                newMin = tsscores.getMinScore(id)
                if (newMin < minScore):
                    minScore = newMin

                nCnsvSites += tsscores.getCnsvNumSites(id)
                cnsvScoreSum += tsscores.getCnsvScoreSum(id) 
                newMin = tsscores.getCnsvMinScore(id)
                if (newMin < cnsvMinScore):
                    cnsvMinScore = newMin

                nNoncnsvSites += tsscores.getNoncnsvNumSites(id)
                noncnsvScoreSum += tsscores.getNoncnsvScoreSum(id) 
                newMin = tsscores.getNoncnsvMinScore(id)
                if (newMin < noncnsvMinScore):
                    noncnsvMinScore = newMin

                    
            if (nSites > 0):
                
                lock.acquire()
                mirna_sites_outfile.write("%s,%d,%3.2f,%3.2f,%d,%3.2f,%3.2f,%d,%3.2f,%3.2f\n" % \
                                              (transcriptID, nSites, scoreSum, minScore, nCnsvSites, cnsvScoreSum, cnsvMinScore,\
                                                   nNoncnsvSites, noncnsvScoreSum, noncnsvMinScore))
                lock.release() 
                
                nSites /= float(seqlen)
                scoreSum /= seqlen
                nCnsvSites /= float(seqlen)
                cnsvScoreSum /= seqlen
                nNoncnsvSites /= float(seqlen)
                noncnsvScoreSum /= seqlen

                lock.acquire() 
                mirna_density_outfile.write("%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n" % \
                                              (transcriptID, nSites, scoreSum, nCnsvSites, cnsvScoreSum, nNoncnsvSites, noncnsvScoreSum))

                lock.release()


    if (args.au_elements):
        # calculate: 
        # 1) count AUUUA
        # 2) number of AU-stretches
        # 3) percentage of UTR covered by AREs 
        # 4) longest A/U stretch 
        auPentamerCount = sequence.count("AUUUA") 

        aurich_re = re.compile(r'[AU]{5,}')  # RE to find more than 5 a/u in a row 
        au_elements = aurich_re.findall(sequence)

        if (au_elements):
            numAUElements = len(au_elements)
            aurichFraction = len("".join(au_elements))/float(seqlen)
            longestAUElement = max(map(len,au_elements))

        else: 
            numAUElements = aurichFraction = longestAUElement = 0

        lock.acquire()
        au_elements_outfile.write("%s,%s,%s,%3.2f,%s\n" % (transcriptID, auPentamerCount, numAUElements, aurichFraction, longestAUElement))
        lock.release()
        


'''MAIN'''
'''if we've been called from the command line and we're in main'''
if __name__ == "__main__":

    '''create a new thread lock'''
    lock = threading.Lock()

    processed = 0 

    delay = 0.1
    if (args.lfold):  # different sleep delays for lfold or regular because threads have different average durations in the two cases
        delay = 0.01

    for line in stdout_from_command(cmd):

        '''decode line'''
        line = line.decode("utf-8")


        if (threading.active_count() < args.nt + 1):
            t = threading.Thread(target=process_line, args=(line,lock,))
            t.start()
            processed += 1 
            
        else:
            while(threading.active_count() >= args.nt + 1):
                time.sleep(delay)
            t = threading.Thread(target=process_line, args=(line,lock,))
            t.start()
            processed += 1 

        if (not(processed % 2500)):
            print ("Processed %d entries..." % processed )

    print ("Processed %d entries." % processed )

    has_stderr = False

    with open ("%s_stderr.log" % args.output, "r") as logfile:
        for i,l in enumerate(logfile):
            has_stderr = True
    if (has_stderr):
        print ("WARNING: stderr was generated during this run. See %s_stderr.log for details" % args.output)
