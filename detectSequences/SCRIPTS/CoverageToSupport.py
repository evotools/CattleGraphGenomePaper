def get_support_value(annot, groups):
    annot = [i for i in annot.split("#") if i in groups]
    values = dict([ i.split("=") for i in annot ])

    return 0

def main():
    import sys

    groups = sys.argv[2].split(",")


    bps = []
    chroms = []
    supp_vecs = []
    # supp_val = []
    repetitive = []
    gap = []
    strand = []
    rlen = []

    if sys.argv[1] == "stdin" or sys.argv[1] == "-":
        dt = sys.stdin
    else:
        dt = open(sys.argv[1])

    if sys.argv[3] == "stdout":
        dr = sys.stdout
    else:
        dr = open(sys.argv[3], "w")

    # Screen for duplicate self alignments withing chromosome
    for line in dt:
        line = line.strip().split()
        chrom, bpi, annot = line[0], line[1], line[3]
        chroms += [chrom]
        bps += [bpi]
        supp_vecs += [ "".join(["1" if i in annot else "0" for i in groups]) ]
        # supp_val += [ ]
        repetitive += [ "1" if "ismask=True" in annot else "0"]
        gap += [ "1" if "isN=True" in annot else "0" ]
        strand += [ "+" if "strand=+" in annot else "-" ]
        rlen += [int(i.split("#")[-1].split("=")[1]) for i in annot.split("#") if "LEN=" in i]

    for i, chrom in enumerate(chroms):
        bp = bps[i] if strand[i] == "+" else rlen[i] - 1 - int(bps[i])
        supp_vec = supp_vecs[i]
        isrep = repetitive[i]
        isgap = gap[i]
        dr.write( "{0}\t{1}\t{1}\t{2}\t{3}\t{4}\n".format(chrom, bp, supp_vec, isrep, isgap) )
    return 0



if __name__ == "__main__":
    main()
    pass