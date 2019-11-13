import sys

if __name__ == "__main__":
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

    for line in ropn:
        if '.gz' in invcf: line = line.decode()
        if "##" in line:
            outfile.write(line)
        elif '##' not in line and '#' in line:
            outfile.write(line)
            ninds = len(line.strip().split()) - 9
        else:
            line = line.strip().split()
            if line[2] != '.': line[2] = '.'
            if line[5] == '.': line[5] = '99'
            if line[6] != '.': line[6] = '.'
            if line[9] != 'GT':
                fmt = line[8].split(':')
                PGTIDX = GTIDX = line[8].split(":").index("GT")
                if 'PGT' in line[8]:
                    PGTIDX = line[8].split(":").index("PGT")
                idGT = [
                    l.split(":")[PGTIDX]
                    if len(l.split(":"))==len(line[8].split(":")) and "|" in l.split(":")[PGTIDX] else l.split(":")[GTIDX]
                    for l in line[9:]
                ]
                line[9:] = idGT
            
            info = line[7].split(';')

            newNfo = {'NA': str(len(line[4].split(','))), 
                        'NS': str(sum([1 for gt in idGT if gt != './.' and gt != '.|.'])),
                        'LEN': '',
                        'TYPE': '',
                        'AC': str(sum([gt.count('1') for gt in idGT]))
                    }
            for f in info:
                if 'SVLEN' in f: newNfo['LEN'] = f.split('=')[1]
                if 'SVTYPE' in f: newNfo['TYPE'] = f.split('=')[1]

            if newNfo['LEN'] == '' and len(line[3]) == len(line[4]) and len(line[3]) == 1:
                newNfo['LEN'] = str(len(line[3]))
                newNfo['TYPE'] = 'snp'
            elif newNfo['LEN'] == '' and len(line[3]) == len(line[4]) and len(line[3]) > 1:
                newNfo['LEN'] = str(len(line[3]))
                newNfo['TYPE'] = 'mnp'
            elif newNfo['LEN'] == '' and len(line[3]) != len(line[4]) and len(line[3]) > len(line[4]):
                newNfo['LEN'] = str(abs(len(line[3]) - len(line[4])))
                newNfo['TYPE'] = 'del'
            elif newNfo['LEN'] == '' and len(line[3]) != len(line[4]) and len(line[3]) < len(line[4]):
                newNfo['LEN'] = str(abs(len(line[3]) - len(line[4])))
                newNfo['TYPE'] = 'ins'

            line[7] = ';'.join(['{}={}'.format(f, newNfo[f].lower()) for f in ['AC', 'LEN', 'NA', 'NS', 'TYPE']])
            outfile.write('{}\n'.format('\t'.join(line)))
    sys.stderr.write('All done.\n\n')

    





