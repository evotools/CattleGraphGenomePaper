



def main():
    import sys
    import os
    if ".gz" in sys.argv[1]:
        import gzip
        opn = gzip.open
        inf = opn(sys.argv[1])
    elif sys.argv[1] == "-":
        inf = sys.stdin.readlines()
    else:
        opn = open
        inf = opn(sys.argv[1])

    ABs = {}
    count = {} 
    QUALs = {}
    DPs = {}

    for n, line in enumerate(inf):
        if n == 0: continue
        CHROM, POS, REF, ALT, QUAL, zQUAL, GT, LEN, GROUP, DP, AD, AB = line.strip().split()
        if GROUP not in ABs:
            ABs[GROUP] = float(AB)
            QUALs[GROUP] = float(zQUAL)
            DPs[GROUP] = float(DP)
            count[GROUP] = 1
            continue
        ABs[GROUP] += float(AB)
        QUALs[GROUP] += float(zQUAL)
        DPs[GROUP] += float(DP)
        count[GROUP] += 1

    Groups = ABs.keys() 
    Groups_sort  = [i for i in Groups if "<" in i] 
    Groups_sort += list(map(str, sorted( [int(i) for i in Groups if ">" not in i and "<" not in i] )))
    Groups_sort += [i for i in Groups if ">" in i]
    print("GROUP,avgAB,avgQUAL,avgDP,N")
    for g in Groups_sort:
        print("{},{},{},{},{}".format(g, ABs[g] / count[g], QUALs[g] / count[g], DPs[g] / count[g], count[g]))

    return 0


if __name__ == "__main__":
    main()


