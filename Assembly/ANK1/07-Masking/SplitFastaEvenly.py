import sys
import os

fname = sys.argv[1]
maxlen = int(sys.argv[2])

if not os.path.exists("./LIST"): os.mkdir("./LIST")

n = 0
c = 1
lines = []
fnames = open("splitList.txt", "w")
lname = open("splitSeqLengths.txt", "w")
for line in open(fname):
        length = int(line.strip().split()[1])
        lines.append(line)
        n += length
        if n > maxlen:
                print >>fnames, "./LIST/ctgList_{}.txt".format(c)
                print >>lname, "./LIST/ctgList_{}.txt\t{}".format(c, n)
                of = open("./LIST/ctgList_{}.txt".format(c), "w")
                for l in lines:
                        print >>of, "{}".format(l.strip().split()[0])
                of.close()
                c += 1
                n = 0
                lines = []

print >>fnames, "./LIST/ctgList_{}.txt".format(c)
of = open("./LIST/ctgList_{}.txt".format(c), "w")
for l in lines: print >>of, "{}".format(l.strip().split()[0])
of.close()
c += 1
n = 0
lines = []
fnames.close()
print "Done"   