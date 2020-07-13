import sys



def parseall():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indata", metavar = 'input.txt.gz', type = str, help = 'Input file',\
                             dest = 'indata', required = True)
    return parser.parse_args()

zval = lambda x,m,s: (x-m)/s

def main():
    import sys, gzip, numpy 
    from scipy import stats
    args = parseall()

    if ".gz" in args.indata:
        opn = gzip.open(args.indata)
    elif args.indata == "-":
        opn = sys.stdin
    else:
        opn = open(args.indata)
    
    QUALS = []
    for n,line in enumerate(opn):
        if n > 0: 
            CHROM, POS, REF, ALT, QUAL, GT, LEN, GROUP, DP, AD, AB = line.strip().split()
            QUALS.append(float(QUAL))
    zQUALS = stats.zscore(numpy.array(QUALS))
    opn.close()

    if ".gz" in args.indata:
        opn = gzip.open(args.indata)
    elif args.indata == "-":
        opn = sys.stdin
    else:
        opn = open(args.indata)
    for n,line in enumerate(opn):
        if n == 0:
            print("CHROM\tPOS\tREF\tALT\tQUAL\tzQUAL\tGT\tLEN\tGROUP\tDP\tAD\tAB")
        if n > 0: 
            CHROM, POS, REF, ALT, QUAL, GT, LEN, GROUP, DP, AD, AB = line.strip().split()
            zQUAL = zQUALS[n-1]
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(CHROM, POS, REF, ALT, QUAL, zQUAL, GT, LEN, GROUP, DP, AD, AB))
    




if __name__ == "__main__":
    main()
