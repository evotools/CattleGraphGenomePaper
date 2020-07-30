import sys


def main():
    for n,line in enumerate(sys.stdin):
        if n == 0: continue
        CHROM, POS, REF, ID = line.strip().split()
        print( "{}\t{}\t{}\t{}".format(CHROM, int(POS) - 1, int(POS) - 1 + len(REF), ID) )
        
if __name__ == "__main__":
    main()