

def parseall():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indata", metavar = 'input.txt.gz', type = str, help = 'Input file',\
                             dest = 'indata', required = True)
    parser.add_argument("-b", "--boundaries", metavar = 'y/N', type = str, help = 'Define boundaries of lenght distibution.',\
                             dest = 'bnd', required = False, default = "-60,60")
    return parser.parse_args()



def isint(val):
    try:
        int(val)
        return True
    except:
        return False


def parse_bnd(bnd):
    bnd = bnd.split(",")
    if len(bnd) == 2:
        if isint(bnd[0]):
            lowerB = int(bnd[0])    
        else:
            lowerB = -999999999
        if isint(bnd[1]): 
            upperB = int(bnd[1])    
        else:
            upperB = 999999999
        return lowerB, upperB
    elif len(bnd) == 1:
        if isint(bnd[0]):
            b = int(bnd[0])
        else:
            return -999999999, 999999999 
        if b < 0: 
            return b, 999999999
        else:
            return -999999999, b
    else:
        return -999999999, 999999999




def main():
    import sys, gzip
    args = parseall()
    lowerB, upperB = parse_bnd(args.bnd)

    if ".gz" in args.indata:
        opn = gzip.open(args.indata)
    elif args.indata == "-":
        opn = sys.stdin
    else:
        opn = open(args.indata)
    
    print("CHROM\tPOS\tREF\tALT\tQUAL\tGT\tLEN\tGROUP\tDP\tAD\tAB")
    for line in opn:
        if line[0] == "#": continue  
        CHR, BP, REF, ALTS, QUAL, GT, DP, AllCounts = line.strip().split()
        geno = GT.split("/")
        AllCount = AllCounts.split(",")
        geno = map(int, geno)
        AllCount = map(float, AllCount)
        A1 = geno[0]
        A2 = geno[1]
        LEN = [ len(i)-len(REF) for i in ALTS.split(",") ]
        GROUP = 0
        for L in LEN:
            if abs(L) > abs(GROUP): GROUP = L
        if GROUP > upperB: GROUP = ">{}".format(upperB)
        elif GROUP < lowerB: GROUP = "<{}".format(lowerB)
        if AllCount[geno[0]] == 0: 
            AB = 1
        elif AllCount[geno[1]] == 0:
            AB = 1
        else:
            AB = AllCount[A1] / (AllCount[A1] + AllCount[A2])
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(CHR, BP, REF, ALTS, QUAL, GT, ','.join(map(str, LEN)), GROUP, DP, AllCounts, AB))

if __name__ == "__main__":
    main()

