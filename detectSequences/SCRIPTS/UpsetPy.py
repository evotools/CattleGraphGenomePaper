



class Counts:
    def __init__(self, log10 = True):
        self.counts = {}
        self.repetitiveness = {}
        self.groups = []
        self.series = None
        self.aslog10 = log10


    def getGroups(self, groups):
        import itertools
        self.groups = groups.split(",")
        lst = list(itertools.product([0, 1], repeat=len(self.groups)))
        self.counts.update( {''.join(map(str, it)): 0 for it in lst} )
        self.repetitiveness.update( {''.join(map(str, it)): 0 for it in lst} )

    def readBed(self, bed, limit = None):
        for line in open(bed):
            chrom, bpi, bpe, annot = line.strip().split()
            try:
                supVec, isRepetitive, nBases = annot.split("#")
            except:
                supVec, isRepetitive, isNmer, nBases = annot.split("#")
            if limit is None:
                self.counts[supVec] += int(bpe) - int(bpi) + 1
            elif limit == isRepetitive:
                self.counts[supVec] += int(bpe) - int(bpi) + 1
            if int(isRepetitive) == "1":
                self.repetitiveness[supVec] += int(nBases)

    def readCov(self, cov, limit = None):
        for line in open(cov):
            try:
                chrom, bp, supVec, isRepetitive = line.strip().split()
            except:
                chrom, bp, supVec, isRepetitive, = line.strip().split()
            if limit is None or limit == isRepetitive:
                self.counts[supVec] += 1
            if int(isRepetitive) == "1":
                self.repetitiveness[supVec] += 1

    def returnTable(self):
        from upsetplot import from_memberships
        vectors, values = [], []
        for vector in self.counts:
            if self.aslog10:
                from math import log
                try: 
                    newvalue = log(self.counts[vector], 10)
                except:
                    newvalue = 1
                values.append( newvalue )
            else:
                values.append( self.counts[vector] )
            classes = [ self.groups[n] for n,i in enumerate(vector) if i == "1" ]
            print(vector, classes, self.counts[vector])
            vectors += [ classes ]
        self.series = from_memberships(vectors, data = values)
        return self.series 

    def saveTable(self, outname, myTable = None):
        if self.series is None and myTable is None:
            self.returnTable()
        if myTable is not None:
            myTable.to_csv("{}.csv".format(outname), header=True)
        elif self.series is not None:
            self.series.to_csv("{}.csv".format(outname), header=True)


    def getPlot(self, outname):
        from upsetplot import plot
        from matplotlib import pyplot
        if self.series is None:
            self.returnTable()
        plot(self.series)
        pyplot.savefig('{}.pdf'.format(outname), dpi = 300)



class Main:

    def _init__(self):
        self.version="0.1.0"
        self.arguments = None
        self.fileType = "guess"
        self.myTable = None
        

    def args(self):
        import argparse
        # Argument definition
        parser = argparse.ArgumentParser()
        parser.add_argument("-i", "--input", metavar = 'coverage.cov/bed', type = str, help = 'Input coverage regions/bases',\
                                dest = 'infile', required = True)
        parser.add_argument("-g", "--groups", metavar = 'class1,class2,...,classN', type = str, help = 'Ordered, comma-separated list of groups',\
                                dest = 'mygroups', required = True)
        parser.add_argument("-t", "--type", metavar = 'bed/cov', type = str, help = 'Input file type, if not specified will guess',\
                                dest = 'ftype', required = False, default = "guess", choices = ["bed", "cov", "guess"])
        parser.add_argument("-o", "--out", metavar = 'prefix', type = str, help = 'Output file prefix.',\
                                dest = 'outname', required = False, default = "upsetplot")  
        parser.add_argument("-l", '--log10', help = 'Show values as log10 (good for large genomes)',\
                                action='store_true')   
        parser.add_argument("-r", '--repetitive', help = 'Analyse only repetitive elements',\
                                action='store_true')   
        parser.add_argument("-n", '--nonrepetitive', help = 'Analyse only non-repetitive elements',\
                                action='store_true')   
        parser.add_argument("-c", '--countNs', help = 'Include Ns in the creation of intervals (not working yet)',\
                                action='store_true')  
        self.arguments = parser.parse_args()


    def guessInputType(self, inputname):
        opn = open
        if ".gz" in inputname:
            import gzip
            opn = gzip.open
        line = opn(inputname).readline()
        line = line.strip().split()
        if len(line) > 4 or len(line) < 4:
            import sys
            sys.exit("File is not compliant.\nPlease check your input files and try again.")
        if "#" in line[3] and line[1].isdigit() and line[2].isdigit():
            self.fileType = "bed"
        elif line[1].isdigit() and line[2].isdigit() and line[3].isdigit():
            self.filetype = "cov"
        else:
            import sys
            sys.exit("File is not compliant.\nPlease check your input files and try again.")


    def checkModules(self):
        installed = [0,0,0]
        modules = ["itertools", "matplotlib", "upsetplot"]
        try:
            import itertools
            installed[0] = True
        except:
            print("No itertools installed")
        try:
            import matplotlib
            installed[1] = True
        except:
            print("No matplotlib installed")
        try:
            import upsetplot
            installed[2] = True
        except:
            print("No upsetplot installed")
        if sum(installed) < len(installed):
            import sys
            joined = ','.join([ i for n,i in enumerate(modules) if installed[n] == 0 ])
            sys.exit("Install {} and try again".format( joined ) )


    def run(self):
        self.checkModules()
        self.args()

        limitTo=None
        if self.arguments.repetitive:
            limitTo="1"
        if self.arguments.nonrepetitive:
            limitTo="0"

        mycounts = Counts(log10 = self.arguments.log10)
        mycounts.getGroups(self.arguments.mygroups)
        if self.arguments.ftype == "guess":
            guessInputType(self.arguments.infile)
        else:
            self.fileType = self.arguments.ftype
        if self.fileType == "bed":
            mycounts.readBed(self.arguments.infile, limit=limitTo)
        else:
            mycounts.readCov(self.arguments.infile, limit=limitTo)
        mycounts.getPlot(self.arguments.outname)
        mycounts.saveTable(self.arguments.outname)
        print("All Done")




if __name__ == "__main__":
    main = Main()
    main.run()
    pass

