import sys


def main():
    for line in sys.stdin:
        if "#" in line[0]: continue
        line = line.strip().split()
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, samples = line.strip().split(8)
        
        print( "{}\t{}\t{}".format(CHROM, POS, END) )


if __name__ == "__main__":
    main()