import sys
import numpy as np

def processSample(bedfile, karyotype, positions):
    # Process single sample
    for line in open(bedfile):
        chromosome, bpi, bpe = line.strip().split()[0:3]
        index = np.where( ( positions[chromosome] >= int(bpi) ) & ( positions[chromosome] <= int(bpe) ) )
        karyotype[chromosome][index] += 1.0

    return karyotype


def generateKaryotype(positions):
    karyotype = {}
    for chromosome in positions:
        karyotype[chromosome] = np.zeros(positions[chromosome].shape[0])
    return karyotype

def getROHclass(length):
    classes = [100000, 1000000]
    if length < classes[0]:
        return "1"
    elif length >= classes[0] and length < classes[1]:
        return "2"
    if length >= classes[1]:
        return "3"

def parsing():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", metavar = 'input.fa', type = str, help = 'Fasta file to process',\
                             dest = 'fasta', required = True)
    parser.add_argument("-o", "--output", metavar = 'outname', type = str, help = 'Output file name',\
                             dest = 'out', required = False, default = "consensus")
    return parser.parse_args()



def seqHighIntervals(listOfBool):
    import re
    intervals = []
    matches = list(re.finditer('1+', listOfBool))
    intervals = [[match.start(), match.end()] for match in matches]
    return intervals

def getConsensus(seqID, seq):
    topscoring = ''.join([ "1" if i.islower() else "0" for i in list(seq) ])
    intervals_temporary = seqHighIntervals(topscoring)
    return [ [seqID] + itv + [ "LEN={}".format(itv[1] - itv[0]) ] for itv in seqHighIntervals(topscoring) ]

def readFasta(infa):
    seqs = {}
    tmpseq = None
    tmpchr = None
    for l in infa:
        if ">" in l:
            if tmpchr is not None: 
                seqs[tmpchr] = tmpseq
            tmpchr = l.replace(">", "").strip().split()[0]
            tmpseq = ""
            continue
        tmpseq += l.strip()
    seqs[tmpchr] = tmpseq
    return seqs
    

def readFastaGz(infa):
    seqs = {}
    tmpseq = None
    tmpchr = None
    for l in infa:
        l = l.decode("utf-8") 
        if ">" in l:
            if tmpchr is not None: 
                seqs[tmpchr] = tmpseq
            tmpchr = l.replace(">", "").strip().split()[0]
            tmpseq = ""
            continue
        tmpseq += l.strip()
    seqs[tmpchr] = tmpseq
    return seqs

if __name__ == "__main__":

    args = parsing()
    # Read positions
    if ".gz" in args.fasta:
        import gzip as gz
        fasta_file = gz.open(args.fasta)
        sys.stderr.write("Reading fasta.gz...\n")
        sequences = readFastaGz(fasta_file)
        sys.stderr.write("Read fasta\n")
    else:
        fasta_file = open(args.fasta)
        sys.stderr.write("Reading fasta...\n")
        sequences = readFasta(fasta_file)
        sys.stderr.write("Read fasta\n")
    
    # Compute ratio and define highly homozygote regions
    sys.stderr.write("Generating masked regions...\n")
    regions = []
    for sequence in sequences:
        regions += getConsensus(sequence, sequences[sequence])

    # Save consensus homozygote regions
    of = open("{}".format(args.out), "w")
    for region in regions:
        of.write("{}\t{}\t{}\t{}\n".format(region[0], region[1], region[2], region[3]))
    sys.stderr.write("Output printed\n")
    of.close()
