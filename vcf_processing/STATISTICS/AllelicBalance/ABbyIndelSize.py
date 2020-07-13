fAB = lambda x,y: x/(x+y)

def main():
    import sys, gzip
    if ".gz" in sys.argv[1]:
        opn = gzip.open(sys.argv[1])
    elif sys.argv[1] == "-":
        opn = sys.stdin
    else:
        opn = open(sys.argv[1])
    print("CHROM\tPOS\tREF\tALT\tGT\tLEN\tAD\tAB")
    for line in opn:
        if line[0] == "#": continue  
        CHR, BP, REF, ALTS, GT, AllDepths = line.strip().split()
        geno = list(map( int, GT.replace("|", "/").split("/") ))
        ADs = AllDepths.split(",")
        ADs = [  float( ADs[i] ) for i in geno ]
        LEN = max( [ len(i) for i in ALTS.split(",") ] + [ len(REF) ] ) - 1
        if len(REF) > max( [ len(i) for i in ALTS.split(",") ] ):
            LEN = LEN * -1
        AB = fAB(ADs[0], ADs[1])
        try:
            AB = fAB(ADs[0], ADs[1])
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(CHR, BP, REF, ALTS, GT, LEN, AllDepths, AB))
        except:
            print(line)
            sys.exit()

if __name__ == "__main__":
    main()
