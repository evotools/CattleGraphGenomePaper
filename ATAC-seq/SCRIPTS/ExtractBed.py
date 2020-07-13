import sys

def newRange(positions, bpiR, bpeR):
	Sbpi, Sbpe = list(map(int, positions.split("-") ))
	bpi = Sbpi + int(bpiR)
	bpe = bpi + (int(bpeR) - int(bpiR))
	return bpi, bpe

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

	outBroad = open("{}.bed".format(outPrefix), "w")
	outNarrow = open("{}.narrowRegion.bed".format(outPrefix), "w")
	regions = {}
	print("Import regions and save narrow peaks...")
	for line in infile:
		line = line.strip().split()
		peakID = line[3]
		if line[0][0] == ":":
			positions = line[0].split(":")[-1]
			sequenceID = ':'.join( line[0].split(":")[2:-1] )		
			Sbpi, Sbpe = list(map(int, positions.split("-") ))
			try: regions["{}#:_:#{}#:_:#{}".format(sequenceID, Sbpi, Sbpe)] += [peakID]
			except: regions["{}#:_:#{}#:_:#{}".format(sequenceID, Sbpi, Sbpe)] = [peakID]
			Nbpi, Nbpe = newRange(positions, line[1], line[2])
			outNarrow.write("{}\t{}\t{}\t{}\n".format(sequenceID, Nbpi, Nbpe, peakID))
		else:
			sequenceID = line[0]
			try: regions["hereford.{}#:_:#{}#:_:#{}".format(sequenceID, line[1], line[2])] += [peakID]
			except: regions["hereford.{}#:_:#{}#:_:#{}".format(sequenceID, line[1], line[2])] = [peakID]
			outNarrow.write("hereford.{}\t{}\t{}\t{}\n".format(sequenceID, line[1], line[2], peakID))	

	# Save broad regions
	print("Save broad regions...")
	for rname in sorted([ k  for k in regions.keys() ]):
		sequence, bpi, bpe = rname.split('#:_:#')
		outBroad.write( "{}\t{}\t{}\t{}\n".format(sequence, bpi, bpe, '#'.join(regions[rname])), ) 

	outBroad.close()
	outNarrow.close()
	print("All done.")

if __name__ == "__main__":
	main()
