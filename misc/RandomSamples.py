

def main():
    import sys, random
    invcf = sys.argv[1]
    otvcf = sys.argv[2]
    ropn = open
    wopn = open
    if '.gz' in invcf:
        import gzip
        ropn = gzip.open(invcf)
    elif "-" == invcf or "stdin" == invcf:
        ropn = sys.stdin
    else:
        ropn = open(invcf)
    if '.gz' in otvcf:
        sys.stderr.write('Output can be only a plain vcf.\nOptionally, you can pipe to bgzip using:\n')
        sys.stderr.write('python GraphVCF.py invcf.vcf(.gz) stdout | bgzip -c > graphvcf.vcf\nSaving as regular vcf...\n')
        otvcf = otvcf.replace('.gz', '')
    elif otvcf.lower() == 'stdout':
        outfile = sys.stdout
    else:     
        outfile = wopn(otvcf, 'w')

    nsamples = int(sys.argv[3])

    SIDs = ["Sample_{}".format(n) for n in range(0, nsamples)]

    for line in ropn:
        if '.gz' in invcf: line = line.decode()
        if "##" in line:
            outfile.write(line)
        elif '##' not in line and '#' in line:
            line = line.strip().split()
            if "FORMAT" not in line:
                line.append("FORMAT")

            line = '\t'.join(line[0:9] + SIDs)
            outfile.write(line + '\n')
        else:
            line = line.strip().split()
            if len(line) >= 9:
                line[8] = "GT"
            else:
                line.append("GT")
            nalleles = 1 + len( line[4].split(",") )
            combinations = ["{}|{}".format(i,k) for i in range(0, nalleles) for k in range(0, nalleles)]
            genotypes = combinations + [random.sample(combinations, 1)[0] for i in range(0,nsamples - len(combinations))]
            line = '\t'.join(line[0:9] + genotypes)
            outfile.write(line + '\n')
    return 0



if __name__ == "__main__":
    main()