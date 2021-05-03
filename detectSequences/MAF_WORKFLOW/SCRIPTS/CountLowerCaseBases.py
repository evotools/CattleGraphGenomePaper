import sys

counts = [0, 0]
tot = [0, 0]
seqID = ''
print( "SEQID\tBPmasked\tBPtot\tPercMasked" )
for line in open(sys.argv[1]):
    if ">" in line:
        if seqID != '':
            print( "{}\t{}\t{}\t{}".format(seqID, counts[1], counts[0], float(counts[1]) / float(counts[0])) )
            tot[0] += counts[0]
            tot[1] += counts[1]
        seqID = line.strip().replace(">", "").split()[0]
        counts = [0, 0]
        continue
    seq = line.strip()
    counts[0] += len(seq)
    counts[1] += sum(map(str.islower, seq))
print( "{}\t{}\t{}\t{}".format("Total", tot[1], tot[0], float(tot[1]) / float(tot[0])) )