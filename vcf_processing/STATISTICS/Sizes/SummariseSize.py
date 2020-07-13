
def main():
    import sys

    classes = {(0, 1): [0, 0, 0, 0, 0],
    (1, 31): [0, 0, 0, 0, 0],
    (31, 101): [0, 0, 0, 0, 0],
    (101, 501): [0, 0, 0, 0, 0],
    (501,1001): [0, 0, 0, 0, 0],
    (1001, 5001): [0, 0, 0, 0, 0],
    (5001,100000001):[0, 0, 0, 0, 0]}

    vects = [0, 1, 31, 101, 501, 1001, 5001, 100000001]
    for line in sys.stdin.readlines():
        if "CHROM" in line: continue
        CHROM, POS, REF, ALT, QUAL, zQUAL, GT, LEN, GROUP, DP, AD, AB = line.strip().split()
        LEN = max([abs(int(i)) for i in LEN.split(",") ])
        GT = [int(i) for i in GT.replace("|","/").split("/")]
        cl = [(vects[i], vects[i + 1]) for i in range(0, len(vects) - 1) if vects[i]<=LEN and vects[i + 1]>LEN ][0]
        classes[ cl ][0] += 1
        if GT[0] == GT[1] and 0 in GT:
            classes[cl][1] += 1
        if GT[0] != GT[1] and 0 in GT:
            classes[cl][2] += 1
        if GT[0] == GT[1] and 0 not in GT:
            classes[cl][3] += 1
        if GT[0] == GT[1] and 0 not in GT:
            classes[cl][4] += 1
        
    print("SIZE,N_VAR,N_RR,N_RA,N_AA_HOM,N_AA_HET")
    for i in range(0, len(vects) - 1):
        cl = (vects[i], vects[i + 1])
        vals = [str(i) for i in classes[cl]]
        if cl[1] > 1000000: 
            print( "{},{}".format(">={}".format(cl[0]), ','.join(vals)) )
        else:
            print( "{},{}".format(cl[1] - 1, ','.join(vals)) )

if __name__ == "__main__":
    main()
