

def main():
    import sys
    infile = sys.argv[1]
    if ".gz" in infile:
        import gzip as gz
        opn = gz.open(infile)
    elif infile == "-":
        opn = sys.stdin
    else:
        opn = open(infile)

    try:
        filt = int(sys.argv[2])
    except:
        filt = None
    
    line1=0
    line2=0
    line3=0

    for n,line in enumerate(opn):
        if "##" in line: 
            if "##FORMAT=<ID=SVLEN" in line: line1=1
            if "##FORMAT=<ID=END" in line: line2=1
            if "##FORMAT=<ID=CIEND" in line: line3=1
            print(line.strip())
            continue
        elif "#" in line and "##" not in line:
            if not line1: print('''##FORMAT=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">''')
            if not line2: print('''##FORMAT=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">''')
            if not line3: print('''##FORMAT=<ID=CIEND,Number=2,Type=Integer,Description="Breakpoint uncertainty for END">''')
            print(line.strip())
            continue
        line = line.strip().split()
        POS= int(line[1])
        SVLENs = [str( len(A2) - len(line[3]) ) for A2 in line[4].split(",") ]
        if filt is not None and sum([1 for i in SVLENs if int(i) <= filt and int(i) >= -filt ]) == len(SVLENs): continue
        BPEs = [POS + abs(int(SVLEN)) for SVLEN in SVLENs ]
        if len(BPEs) == 1:
            END=BPEs[0]
            CIEND='0,0'
        else:
            END = int( sum(BPEs) / float(len(BPEs)) )
            CI = sorted( [i-END for i in BPEs] )
            CIEND = ','.join( [ str( CI[0] ), str( CI[-1] ) ] )
        SVLEN=','.join(SVLENs)
        infofield = line[7].split(';')
        presentSVL=0
        presentCI=0
        presentED=0
        for n,f in enumerate(infofield):
            if "SVLEN=" in f: infofield[n] = "SVLEN={}".format( SVLEN ); presentSVL=1
            if "END=" in f: infofield[n] = "END={}".format( END ); presentED=1
            if "CIEND=" in f: infofield[n] = "CIEND={}".format( CIEND ); presentCI=1
        line[7] = ';'.join(infofield)
        if not presentSVL: line[7] = "{};SVLEN={}".format( line[7], SVLEN )
        if not presentED: line[7] = "{};END={}".format( line[7], str(END) )
        if not presentCI: line[7] = "{};CIEND={}".format( line[7], CIEND )
        print('\t'.join(line))
if __name__ == "__main__":
    main()
