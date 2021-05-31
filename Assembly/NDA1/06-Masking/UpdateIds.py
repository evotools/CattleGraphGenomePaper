import sys

codes = {i.strip().split()[0]:i.strip().split()[1] for i in open(sys.argv[2])}

ncodes = open("FastaCodeCorrespondence.txt", "w")
nwname = sys.argv[1].replace(".fasta", "").replace(".fa","").replace(".fna","") + ".renamed.fa"
nfasta = open(nwname, "w")
c = 1
for line in open(sys.argv[1]):
        if ">" in line:
                name = line.strip().replace(">","")
                name = codes.get(name, name)
                if len(name) > 50 or "[" in name or ":" in name:
                        name = "seq{}".format(c)
                        c += 1
                ncodes.write("{},{}\n".format(line.strip().replace(">",""), name))
                nfasta.write(">{}\n".format(name))
                continue
        nfasta.write(line)