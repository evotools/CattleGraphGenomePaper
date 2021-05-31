import sys

conv={str(n+1): i.strip().split()[0] for n,i in enumerate(open(sys.argv[1]))}

for i in open(sys.argv[2]):
        print('\t'.join([conv.get(i.strip().split()[0], i.strip().split()[0]), i.strip().split()[1], i.strip().split()[2], i.strip().split()[3]]))