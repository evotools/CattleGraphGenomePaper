#!/usr/bin/env python

# Gather contigs statistics after gap filling
#

# Read one fasta file
def readFasta(inputfasta):
    opn = open
    if ".gz" in inputfasta:
        import gzip as gz
        opn = gz.open
    infa = opn(inputfasta)
    seqs = {}
    tmpseq = None
    tmpchr = None
    for l in infa:
        try:
            l = l.decode("utf-8") 
        except:
            l = l
        if ">" in l:
            if tmpchr is not None: 
                seqs[tmpchr] = tmpseq
            tmpchr = l.replace(">", "").strip().split()[0]
            tmpseq = ""
            continue
        tmpseq += l.strip()
    seqs[tmpchr] = tmpseq
    return seqs


# Write output statistics
def statsPrinter(sequences, seqtype, outname, outtype = "text"):
    logging.info(" - {} summary statistics".format(seqtype))
    
    try:
        from prettytable import PrettyTable
        myTable = PrettyTable(["x", "Nx", "Lx", "NGx", "LGx", "BPyield", "GC%", "#Gaps", "nN"])
        myTable.align["x"] = "l"
        myTable.padding_width = 1 # One space between column edges and contents (default) 
        for q in map(str, sequences.qtiles):
            myTable.add_row([q, sequences.Ns[q], sequences.Ls[q], sequences.NGs[q], sequences.LGs[q], sequences.bpYield[q], sequences.Gc[q], sequences.nGaps[q], sequences.nNs[q]])
        logging.info("\n" + myTable)
    except:
        logging.info("No prettytable installed.")
        mLength = max([len(str(sequences.Ns[q])) for q in map(str, sequences.qtiles)] + \
                    [len(str(sequences.Ls[q])) for q in map(str, sequences.qtiles)] + \
                    [len(str(sequences.NGs[q])) for q in map(str, sequences.qtiles)] + \
                    [len(str(sequences.LGs[q])) for q in map(str, sequences.qtiles)] + \
                    [len(str(sequences.bpYield[q])) for q in map(str, sequences.qtiles)]) + 2
        logging.info('|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|'.format("x", "Nx", "Lx", "NGx", "LGx", "BPyield", "GC%", "#Gaps", "nN", width = mLength))
        logging.info("-"* mLength * 9)
        for q in map(str, sequences.qtiles):
            logging.info('|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|'.format(q, sequences.Ns[q], sequences.Ls[q], \
                                                                                                            sequences.NGs[q], sequences.LGs[q], \
                                                                                                            sequences.bpYield[q], \
                                                                                                            sequences.Gc[q], sequences.nGaps[q], \
                                                                                                            sequences.nNs[q], width = mLength))


# Write intervals
def writeIntervals(intervals, outname, outtype):
    ofm = "txt"
    delimiter = ' '
    if outtype == "csv":
        ofm = "csv"
        delimiter = ","
    elif outtype == "csv2": 
        ofm = "csv"
        delimiter = ";"
    elif outtype == "tsv": 
        ofm = "tsv"
        delimiter = "\t"
    of = open("{}_itv.{}".format(outname, ofm), "w")
    for interval in intervals:
        of.write("{}\n".format(delimiter.join(map(str, interval))))
    return 0

# Write output statistics
def statsWriter(sequences, seqtype, outname, outtype = "text"):
    logging.info(" - {} summary statistics".format(seqtype))
    ofm = "txt"
    delimiter = ' '
    if outtype == "csv":
        ofm = "csv"
        delimiter = ","
    elif outtype == "csv2": 
        ofm = "csv"
        delimiter = ";"
    elif outtype == "tsv": 
        ofm = "tsv"
        delimiter = "\t"
    
    with open("{}_{}.{}".format(seqtype, outname, ofm), "w") as sm:
        # Create actual output
        if outtype == "text":
            try:
                from prettytable import PrettyTable
                myTable = PrettyTable(["x", "Nx", "Lx", "NGx", "LGx", "BPyield", "GC%", "#Gaps", "nN"])
                myTable.align["x"] = "l"
                myTable.padding_width = 1 # One space between column edges and contents (default) 
                for q in map(str, sequences.qtiles):
                    myTable.add_row([q, sequences.Ns[q], sequences.Ls[q], sequences.NGs[q], sequences.LGs[q], sequences.bpYield[q], sequences.Gc[q], sequences.nGaps[q], sequences.nNs[q]])
                sm.write(myTable.get_string())
                return 0
            except:
                mLength = max([len(str(sequences.Ns[q])) for q in map(str, sequences.qtiles)] + \
                            [len(str(sequences.Ls[q])) for q in map(str, sequences.qtiles)] + \
                            [len(str(sequences.NGs[q])) for q in map(str, sequences.qtiles)] + \
                            [len(str(sequences.LGs[q])) for q in map(str, sequences.qtiles)]) + \
                            [len(str(sequences.bpYield[q])) for q in map(str, sequences.qtiles)] + 2
        
                sm.write('|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|\n'.format("x", "Nx", "Lx", "NGx", "LGx", "BPyield", "GC%", "#Gaps", "nN", width = mLength))
                for q in map(str, sequences.qtiles):
                    sm.write('|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|\n'.format(q, sequences.Ns[q], sequences.Ls[q], \
                                                                                                            sequences.NGs[q], sequences.LGs[q], \
                                                                                                            sequences.bpYield[q], \
                                                                                                            sequences.Gc[q], sequences.nGaps[q], \
                                                                                                            sequences.nNs[q], width = mLength))
        elif outtype == "csv" or outtype == "csv2" or outtype == "tsv":
            sm.write('{}\n'.format( delimiter.join(["x", "Nx", "Lx", "NGx", "LGx", "BPyield", "GC%", "#Gaps", "nN"])))
            for q in map(str, sequences.qtiles):
                    sm.write('{}\n'.format( delimiter.join( map( str, [q, sequences.Ns[q], sequences.Ls[q], \
                                            sequences.NGs[q], sequences.LGs[q], \
                                            sequences.bpYield[q], \
                                            sequences.Gc[q], sequences.nGaps[q], \
                                            sequences.nNs[q]]) ) ) )
        

# Class to create the sequences object.
class Sequences:
    # Sequences basic parameters
    def __init__(self):
        self.qtiles = []
        self.seqs = []
        self.seqIDs = []
        self.NGs = {}
        self.LGs = {}
        self.Ns = {}
        self.Ls = {}
        self.Gc = {}
        self.nGaps = {}
        self.nNs = {}
        self.bpYield = {}
        self.lengths = []
        self.cumulative = []
        self.tNs = 0
        self.tNgaps = 0
        self.totSeq = 0
        self.rlength = 0
        self.elength = 0
        self.GCcounts = []
        self.gapcount = []
        self.cumulativegap = []
        self.ncount = []
        self.cumulativencount = []
        self.cumulativeGC = []
        self.Eratios = []
        self.Rratios = []

    # Importing, sorting and measuring sequences
    def addSeqs(self, seqlist, exp_length):
        logging.info("Parsing and sorting...")
        self.seqs = seqlist
        self.seqs.sort(key = len, reverse = True)
        self.elength = exp_length
        logging.info("Sequence sorted.")

    def addSeqIds(self, seqIDsList):
        self.seqIDs = seqIDsList


    # Gather base metrics for the sequences (length, GC count and %, number of gaps and length)
    def getMetrics(self):
        logging.info("Gathering base metrics...")
        self.lengths = [len(i) for i in self.seqs]
        self.cumulative = [sum(self.lengths[0:i+1]) for i in range(0, len(self.lengths))]
        # General length stats
        self.totSeq = len(self.seqs)
        self.rlength = sum(self.lengths)
        # General GC stats
        self.GCcounts = [i.lower().count("g") + i.lower().count("c") for i in self.seqs]
        self.cumulativeGC = [sum(self.GCcounts[0:i+1]) for i in range(0, len(self.GCcounts))]
        # General gap counts and stats
        self.gapcount = [self.seqGaps(seq) for seq in self.seqs]
        self.tNgaps = sum(self.gapcount)
        self.cumulativegap = [sum(self.gapcount[0:i+1]) for i in range(0, len(self.gapcount))]
        # Get number of N stats
        self.ncount = [seq.count("N") for seq in self.seqs]
        self.tNs = sum(self.ncount)
        self.cumulativencount = [sum(self.ncount[0:i+1]) for i in range(0, len(self.ncount))]
        # Get estimated and real ratio
        self.Eratios = [(value/self.elength) * 100 for value in self.cumulative]
        self.Rratios = [(value/self.rlength) * 100 for value in self.cumulative]
        logging.info("")
        logging.info("There are {} sequences".format(self.totSeq))
        logging.info("Genome expected length: {}".format(self.elength))
        logging.info("Genome real length: {}".format(self.rlength))
        logging.info("GC count: {}".format(sum(self.GCcounts)))
        logging.info("GC %: {}".format( round(   (sum(self.GCcounts) / float(self.rlength)) * 100, 2) ))
        logging.info("N gaps: {}".format(self.tNgaps))
        logging.info("Total Ns in genome: {}".format(self.tNs))
        logging.info("")

    # Function to get N/Lx stats for the sequences provided
    def getNstats(self):
        logging.info("Calculating Nx/Lx/GCx statistics...")
        for N in self.qtiles:
            n = str(N)
            try:
                self.LGs[n]     = sum([1 for i in self.Eratios if i < N])
                self.NGs[n]     = self.lengths[self.LGs[n]]
                self.Ls[n]      = sum([1 for i in self.Rratios if i < N])
                self.Ns[n]      = self.lengths[self.Ls[n]]
                self.Gc[n]      = round((self.cumulativeGC[self.Ls[n]] / self.cumulative[self.Ls[n]]) * 100, 2)
                self.nGaps[n]   = self.cumulativegap[self.Ls[n]] 
                self.nNs[n]     = self.cumulativencount[self.Ls[n]]
                self.bpYield[n] = self.cumulative[self.Ls[n]]
            except:
                logging.error(str(N) + " not found!")
        self.qtiles += [100]
        n = "100"
        self.LGs[n]     = self.totSeq
        self.NGs[n]     = self.lengths[-1]
        self.Ls[n]      = self.totSeq
        self.Ns[n]      = self.lengths[-1]
        self.Gc[n]      = round((sum(self.GCcounts) / sum(self.lengths)) * 100, 2)
        self.nGaps[n]   = self.tNgaps
        self.nNs[n]     = self.tNs
        self.bpYield[n] = self.cumulative[-1]

    # Create intervals quartiles for Lx/Nx statistics
    def defineQtiles(self, initV, endV, intV):
        self.qtiles = list(range(initV, endV, intV))
        
    # Getting gaps stats for N sequences
    def getGaps(self, nSeqs = 0):
        if nSeqs == 0:
            nSeqs = len(self.seqs)
        nGaps = 0
        totNs = 0
        for seq in self.seqs[0:nSeqs]:
            seq = list(seq)
            nNs = 0
            for c in seq:
                if c == "N": 
                    nNs += 1
                    continue
                elif c != "N" and nNs != 0:
                    nGaps += 1
                    totNs += nNs
                    nNs = 0
                    continue
                else:
                    continue
            if nNs != 0:
                nGaps += 1
                totNs += nNs
                nNs = 0
        return totNs, nGaps

    # Getting number of gaps per single sequence
    def seqGapsLegacy(self, seq):
        nGaps = 0
        isGap = False
        seq = list(seq)
        for c in seq:
            if c == "N": 
                isGap = True
                continue
            elif c != "N" and isGap:
                nGaps += 1
                isGap = False
                continue
            else:
                continue
        if isGap:
            nGaps += 1
        return nGaps

    def seqGaps(self, seq):
        import re
        matches = len(list(re.finditer('N+', seq)))
        return matches

    
    def seqGapsIntervals(self, fastaObj):
        import re
        seqIDs = fastaObj.keys()
        if len(seqIDs) == 0:
            return 1
        intervals = []
        for seqID in seqIDs:
            seq = fastaObj[seqID]
            matches = list(re.finditer('N+', seq))
            intervals += [[seqID, match.start() + 1, match.end(), "+", match.end() - match.start() + 1, "Nmer"] for match in matches]
        return intervals

def parse_length(l):
    if l.isdigit():
        l = int(l)
    elif 'g' in l.lower():
        l = int(float(l[0:-1]) * 1000000000)
    elif 'm' in l.lower():
        l = int(float(l[0:-1]) * 1000000)
    elif 'k' in l.lower():
        l = int(float(l[0:-1]) * 1000)
    return(l)


if __name__ == "__main__":
    import sys
    import argparse
    import logging
    import gzip

    # Checking python version
    if sys.version_info[0] < 3:
        logging.error("Python version must be 3.X")
        raise Exception()

    # Argument definition
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", metavar = 'fastafile.fa(.gz)', type = str, help = 'Input fasta file',\
                             dest = 'fasta', required = True)
    parser.add_argument("-l", "--length", metavar = '2700000000', type = str, help = 'Expected genome length',\
                             dest = 'elength', required = True)
    parser.add_argument("-i", "--interval", metavar = '5,100,5', type = str, help = 'Intervals of N/L to consider (default, intervals from 5 to 100 with step = 5)',\
                             dest = 'itv', required = False, default = "5,100,5")
    parser.add_argument("-s", "--split", metavar = 'Y/n', type = str, help = 'Split scaffold in contigs by Ns',\
                             dest = 'do_split', required = False, default = "y")
    parser.add_argument("-o", "--out", metavar = 'prefix', type = str, help = 'Output file prefix.',\
                             dest = 'outname', required = False, default = "stats")
    parser.add_argument("-p", "--positions", metavar = 'y/N', type = str, help = 'Output file with Nmer intervals.',\
                             dest = 'printItvs', required = False, default = "n")
    parser.add_argument("-t", "--type", metavar = 'text/csv/csv2/tsv', type = str, help = 'Output file type (text, comma-sep, semicolumn-sep or tsv).',\
                             dest = 'outtype', required = False, default = "text", choices = ["text", "csv", "csv2", "tsv"])
    args = parser.parse_args()


    # create logger with 'spam_application'
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='{}.log'.format(args.outname),
                    filemode='w')
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    # Print software version
    logging.info("")
    logging.info("ABS: Assembly Base Stats")
    logging.info("")
    logging.info("Read input fasta...")
    seqs = readFasta(args.fasta)
    seqIDs = seqs.keys()
    elength = parse_length(args.elength)
    intervals = list(map(int, args.itv.split(",")))
    

    # Get statistics for scaffolds
    logging.info("Processing scaffolds")
    scaffolds = Sequences()
    scaffolds.addSeqs([seqs[seqID] for seqID in seqIDs], elength)
    scaffolds.addSeqIds(seqIDs)
    scaffolds.defineQtiles(intervals[0], intervals[1], intervals[2])
    scaffolds.getMetrics()
    scaffolds.getNstats()
    if args.printItvs.lower() == "y":
        intervals = scaffolds.seqGapsIntervals(seqs)
        writeIntervals(intervals, args.outname, args.outtype)
    statsPrinter(scaffolds, "scflds", args.outname, args.outtype)
    statsWriter(scaffolds, "scflds", args.outname, args.outtype)
      

    #
    # Start contig level analysis
    #
    newSeq = []
    logging.info("Splitting")
    for seqID in seqs:
        newSeq += [seq for seq in seqs[seqID].split("N") if len(seq) > 0]
    logging.info("Processing contigs")
    intervals = list(map(int, args.itv.split(",")))
    contigs = Sequences()
    contigs.addSeqs(newSeq, elength)
    contigs.defineQtiles(intervals[0], intervals[1], intervals[2])
    contigs.getMetrics()
    contigs.getNstats()
    statsPrinter(contigs, "ctgs", args.outname, args.outtype)
    statsWriter(contigs, "ctgs", args.outname, args.outtype)
    
    if args.do_split.lower() == "y":
        logging.info("Saving split fasta")
        with open("splitFa_{}.fa".format(args.outname), "w") as of:
            n = 1
            for seq in newSeq:
                of.write(">seq{}\n".format(n))
                of.write("{}\n".format(seq))
                n += 1

    #print("Done")
    logging.info("Done")
