import sys

convID=open(sys.argv[1] + ".convid.txt", "w")
n=1
for line in open(sys.argv[1]):
        line = line.strip()
        if ">" in line:
                newID = "seq{}".format(n)
                n += 1
                print(">{}".format(newID))
                convID.write("{}\t{}\n".format(line.replace(">", ""), newID))
        else:
                print(line.upper())