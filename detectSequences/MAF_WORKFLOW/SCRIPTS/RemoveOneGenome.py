
import sys


def processFile(infile, toExclude):
    ofname = open( infile.replace(".bed", ".nonRef.bed"), "w" )
    for line in open(infile):
        myline = line.strip().split()
        if list(myline[-1].split("#")[0])[toExclude] == "0": 
            ofname.write(line)


def main():
    toExclude = int(sys.argv[-1]) - 1
    for infile in sys.argv[1:-1]:
        processFile(infile, toExclude)
        print("Done {}".format(infile))


if __name__ == "__main__":
    main()