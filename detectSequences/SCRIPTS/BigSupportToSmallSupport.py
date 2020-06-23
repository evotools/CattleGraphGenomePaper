

def main():
    import sys
    inputf = sys.argv[1]
    outputbed = sys.argv[2]

    if inputf == "-" or inputf == "stdin":
        indata = sys.stdin
    else:
        indata = open(inputf)

    if outputbed == "stdout":
        bedfile = sys.stdout
    else:
        bedfile = open(outputbed, "w")


    line = indata.readline().strip().split()
    baseSupportVector = "{}#{}#{}".format( line[3], line[4], line[5] ) 
    init_bp = int(line[1])
    end_bp = int(line[2])
    init_chromosome = line[0]
    nbases = 1
    for line in indata:
        chromosome, bp, bpe, supportVector, repetitiveB, gapB = line.strip().split()
        bp = int(bp)
        bpe = int(bpe)
        if "{}#{}#{}".format( supportVector, repetitiveB, gapB ) == baseSupportVector and chromosome == init_chromosome:
            end_bp = bpe
            nbases += 1
        else:
            bedfile.write( "{}\t{}\t{}\t{}#{}\n".format(init_chromosome, init_bp, end_bp+1, baseSupportVector, nbases) )
            nbases = 1
            init_bp, end_bp = bp, bp
            init_chromosome = chromosome
            baseSupportVector = "{}#{}#{}".format( supportVector, repetitiveB, gapB ) 

    bedfile.write( "{}\t{}\t{}\t{}#{}\n".format(init_chromosome, init_bp, end_bp+1, baseSupportVector, nbases) )
    return 0


if __name__ == "__main__":
    main()