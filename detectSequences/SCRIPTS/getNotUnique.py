import sys

for line in open(sys.argv[1]):
    linesplit = line.strip().split()
    supvec = sum(map( int, list(linesplit[-1].split("#")[0]) ))