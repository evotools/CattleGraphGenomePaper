#####
#
# Created by: Andrea Talenti
# Date: 03 Dec, 2019
# 
# Script that create a coverage level from 
# the initial reference-based maf files.
# MAF files have to be created starting from
# the output hal file generated from CACTUS.
# Conversion from hal to maf with hal2maf tool. 
#
#####
import sys
import numpy as np

def getRefQuery(lines, rname, seq, ignoreself):
    if seq is not None:
        references = [line for line in lines if seq in line]
        if ignoreself:
            queries = [line for line in lines if seq not in line and rname not in line]
        else:
            queries = [line for line in lines if seq not in line]
    else:
        references = [line for line in lines if rname in line]
        queries = [line for line in lines if rname not in line]
    return references, queries

def worker(inputs):
    import sys, os, random

    lines, ctr, seq, ignoreSelf = inputs
    seed = random.randint(1000000,9999999)
    fname = "tmp_{}_{}".format(ctr, seed)
    while os.path.exists(fname): fname = "tmp_{}_{}".format(ctr, seed)
    tmpof = open(fname, "w")

    rname = lines[0].strip().split()[1].split(".")[0]
    references, queries = getRefQuery(lines, rname, seq, ignoreSelf)

    for n, reference in enumerate(references):
        newlines = [r for x,r in enumerate(references) if x!=n] + queries
        ref = reference.strip().split()
        refChrom, refBpi, refSeqLength, refseq = ref[1], int(ref[2]), int(ref[3]), np.array(list(ref[6]))
        positions = np.arange(refBpi, refBpi + refSeqLength)
        avail = refseq != "-"
        if sum(avail) != int(refSeqLength):
            refSeqLength = sum(avail)
        refseq = np.array(list(refseq))
        cvgs = {"tot" : np.zeros(refSeqLength).astype("uint8")}
        cvgs.update( {line.split()[1].split(".")[0]: np.zeros(refSeqLength).astype("uint8") for line in lines} )

        for n in range(0,len(newlines)):
            line = newlines[n].strip().split()
            qryId, qryseq = line[1].split(".")[0], np.array(list(line[6]))
            qryseq = qryseq[avail]
            cvgs[qryId] = cvgs[qryId] + ( qryseq != "-" )
            cvgs["tot"] = cvgs["tot"] + ( qryseq != "-" )

        for i in range(0, refSeqLength):
            mystr = "#".join([ "{}={}".format( k, cvgs[k][i] ) for k in cvgs ] )
            tmpof.write( "{0}\t{1}\t{1}\t{2}#ismask={3}#isN={4}\n".format(refChrom, positions[i], mystr, refseq[avail][i].islower(), refseq[avail][i].lower() == "n") )
    if ctr % 10000 == 0:
        sys.stderr.write( "Processed chunk {}\n".format(ctr) )
    tmpof.close()
    return ( ctr,  fname)


def worker_highmem(inputs):
    import numpy as np
    lines, ctr, seq, ignoreSelf = inputs
    rname = lines[0].strip().split()[1].split(".")[0]
    references, queries = getRefQuery(lines, rname, seq, ignoreSelf)
    finalStr = ""
    for n, reference in enumerate(references):
        if not ignoreSelf: newlines = [r for x,r in enumerate(references) if x!=n] + queries
        ref = reference.strip().split()
        refChrom, refBpi, refSeqLength, refseq, refStrand, refTotLen = ref[1], int(ref[2]), int(ref[3]), np.array(list(ref[6])), ref[4], ref[5]
        positions = np.arange(refBpi, refBpi + refSeqLength)
        avail = refseq != "-"
        if sum(avail) != int(refSeqLength):
            sys.exit("Something is wrong with the data.")
        refseq = np.array(list(refseq))
        cvgs = {"tot" : np.zeros(refSeqLength).astype("uint8")}
        cvgs.update( {line.split()[1].split(".")[0]: np.zeros(refSeqLength).astype("uint8") for line in lines} )

        for n in range(0,len(newlines)):
            line = newlines[n].strip().split()
            qryId, qryseq = line[1].split(".")[0], np.array(list(line[6]))
            qryseq = qryseq[avail]
            cvgs[qryId] = cvgs[qryId] + ( qryseq != "-" )
            cvgs["tot"] = cvgs["tot"] + ( qryseq != "-" )

        for i in range(0, refSeqLength):
            mystr = "#".join([ "{}={}".format( k, cvgs[k][i] ) for k in cvgs ] )
            finalStr += "{0}\t{1}\t{1}\t{2}#ismask={3}#isN={4}#strand={5}#LEN={6}\n".format(refChrom, positions[i], mystr, refseq[avail][i].islower(), refseq[avail][i].lower() == "n", refStrand, refTotLen)
    if ctr % 10000 == 0:
        sys.stderr.write( "Processed chunk {}\n".format(ctr) )
    return ( ctr,  finalStr)


def reader(args):
    lines = []
    blocks = []
    infile = open(args.maf)
    for line in infile:
        if "#" in line:
            continue
        elif line[0] == "a": 
            continue
        elif line[0] == "s":
            lines.append(line)
            continue
        elif line.strip() == "" and len(lines) != 0:
            blocks.append(lines)
            lines = []
        else:
            lines = []
    return blocks
    

def spawn_proc(funct, chunk, seq, ignoreself, nthreads):
    if nthreads == 1:
        sys.stderr.write("Start processing the data\n")
        return [funct((block, n, seq, ignoreself)) for n, block in enumerate(chunk)]
    else:
        sys.stderr.write("Start processing the data using {} processes.\n".format(nthreads) )
        from multiprocessing import Pool
        pool = Pool(nthreads)
        return pool.map(funct, zip(chunk,list(range(0, len(chunk))), [seq for i in range(0, len(chunk))], [ignoreself for i in range(0, len(chunk))] ) )


def gather(dat, outfile):
    import os
    for i in dat:
        if os.path.isfile(i[1]):
            [outfile.write(k) for k in open(i[1])]
            os.remove(i[1])
        else:
            outfile.write(i[1])
    return None


def highMemoryPipeline(blocks, args, of):
    if args.seq is not None:
        sys.stderr.write("Limiting to {}\n".format(args.seq))
    if args.ignoreself:
        sys.stderr.write("Ignoring same-genome alignments on different chromosomes\n")
    gather(spawn_proc(worker_highmem, blocks, args.seq, args.ignoreself, args.threads), of)
    return 0


def lowMemoryPipeline(chunks, args, of):
    import os, sys
    if args.seq is not None:
        sys.stderr.write("Limiting to {}\n".format(args.seq))
    if args.ignoreself:
        sys.stderr.write("Ignoring same-genome alignments on different chromosomes\n")
    for w, chunk in enumerate(chunks):
        gather(spawn_proc(worker, chunk, args.seq, args.ignoreself, args.threads), of)
        sys.stderr.write("Processed chunk #{}\n".format(w) )
    return 0


def parse_args():
    import argparse
    # Argument definition
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--maf", metavar = 'maffile.maf', type = str, help = 'Input alignments',\
                            dest = 'maf', required = True)
    parser.add_argument("-t", "--threads", metavar = '1', type = int, help = 'Number of threads to use',\
                            dest = 'threads', required = False)
    parser.add_argument("-s", "--seq", metavar = 'tgtseq', type = str, help = 'Target sequence to compute coverage of',\
                            dest = 'seq', required = False, default = None)
    parser.add_argument("-o", "--out", metavar = 'prefix', type = str, help = 'Output file prefix.',\
                            dest = 'outname', required = False, default = "coverage")
    parser.add_argument("--highmem", help = 'Use high memory instead of memory saving approach',\
                            action='store_true')
    parser.add_argument("--ignoreself", help = 'Ignore self alignments on different chromosomes',\
                            action='store_true')
    
    return parser.parse_args()


def main():
    import os
    import sys
    import time

    args = parse_args()
    # Set environment number of variables
    os.environ["OMP_NUM_THREADS"] = "1"

    sys.stderr.write("Reading alignments blocks and chunking it...\n")
    blocks = reader(args)

    # Create temporary file
    if args.outname == "stdout":
        of = sys.stdout
    else:
        of = open( "{}.bed".format(args.outname), "w" )

    s_time = time.time()
    if not args.highmem:
        sys.stderr.write("Using memory efficient method.\n")
        chunks = [ blocks[i: i + 100000] for i in range(0, len(blocks), 100000) ]
        lowMemoryPipeline(chunks, args, of)
    else:
        sys.stderr.write("Using high memory method.\n")
        highMemoryPipeline(blocks, args, of)

    of.close()

    sys.stderr.write("All done in {}s\n".format(time.time() - s_time))
    return 0


if __name__ == "__main__":
    main()
    