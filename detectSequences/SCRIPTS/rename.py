import sys

of1 = open(sys.argv[2], "w")
conv= open(sys.argv[2] + ".conv", "w")
c=1
for line in open(sys.argv[1]):
    if ">" in line:
        of1.write(">seq{}\n".format(c))
        conv.write("{}\tseq{}\n".format( line.strip().replace(">", ""), c ) )
        c+=1
        continue
    of1.write("{}".format(line))
