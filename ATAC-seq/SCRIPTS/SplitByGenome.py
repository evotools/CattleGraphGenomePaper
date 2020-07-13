import sys


def main():
	filename = sys.argv[1]
	outPrefix = sys.argv[2]
	if filename == "-" or filename == "stdin":
		infile = sys.stdin
	elif ".gz" in filename:
		import gzip as gz
		infile = gz.open(filename)
	else:
		infile = open(filename)

	outdat = {}
	for line in infile:
		line = line.strip().split()
		genome, sequence = line[0].split(".")[0], '.'.join(line[0].split(".")[1:])
		if genome not in outdat:
			outdat[genome] = [ [sequence, line[1], line[2], line[3] ] ]
			continue
		outdat[genome] += [ [sequence, line[1], line[2], line[3] ] ]
	for genome in outdat.keys():
		outfile = open("{}.{}.bed".format( outPrefix, genome), "w")
		for line in outdat[genome]:
			outfile.write("{}\n".format('\t'.join(line)))
		outfile.close()


if __name__ == "__main__":
	main()
