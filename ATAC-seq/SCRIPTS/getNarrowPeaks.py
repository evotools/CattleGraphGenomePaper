import sys

def main():
    narrowPeaks = sys.argv[1]
    bedfile = sys.argv[2]
    peaknames = [line.split()[3] for line in open(bedfile) ]

    for line in open(narrowPeaks):
        if line.strip().split()[3] in peaknames:
            print(line.strip())


if __name__ == "__main__":
    main()