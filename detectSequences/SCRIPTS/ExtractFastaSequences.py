import sys

def main():
    kept = [ line.strip().split()[0] for line in open(sys.argv[2]) ] 
    seqid = ""
    seq_order = []
    sequences = {}
    for line in open(sys.argv[1]):
        if ">" in line:
            line = line.strip().replace(">", "")
            seqid = line
            seq_order += [seqid]
            sequences[seqid] = ""
            continue
        sequences[seqid] += line.strip()

    for s in kept:
        print(">{}".format(s))
        print(sequences[s])

    return 0

if __name__ == "__main__":
    main()