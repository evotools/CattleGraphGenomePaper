import sys

try:
	convID={i.strip().split(",")[1]: i.strip().split(",")[0] for i in open(sys.argv[2]) if "seq" in i.strip().split(",")[1]}
except:
	convID={i.strip().split("\t")[1]: i.strip().split("\t")[0] for i in open(sys.argv[2]) if "seq" in i.strip().split("\t")[1]}
 
for line in open(sys.argv[1]):
	if ">" in line:
		line = line.strip().split(":")[0].replace(">", "")
		line = convID.get(line, line)
		print(">{}".format(line))
		continue
	print(line.strip())