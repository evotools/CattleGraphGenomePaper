import sys


def main():
    buffer = int(sys.argv[1])
    for line in sys.stdin:
        CHROM, POS, END, ID, SVTYPE, SVLEN, CIPOS, CIEND = line.strip().split()
        if int(POS) < int(END):
            CIPOS = min(map(int, CIPOS.split(",")))
            CIEND = max(map(int, CIEND.split(",")))
            POS = int(POS) + CIPOS - buffer - 1
            END = int(END) + CIEND + buffer - 1
        else:
            CIPOS = max(map(int, CIPOS.split(",")))
            CIEND = min(map(int, CIEND.split(",")))
            POS = int(POS) + CIPOS + buffer - 1
            END = int(END) + CIEND - buffer - 1
            POS, END = END, POS
        print( "{}\t{}\t{}\t{}".format(CHROM, POS, END, ID) )


if __name__ == "__main__":
    main()