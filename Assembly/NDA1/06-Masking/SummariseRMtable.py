import sys
import os

def parse_tbl(tblname, dataArray, keyorder):
    myfile = open(tblname)
    start_import = False
    neq = 0
    while not start_import:
        line = myfile.readline()
        if "-----" in line:
            start_import = True
            break
    # Stat actual importation of the table
    for line in myfile:
        if "=" in line: break
        elif len(line.strip()) < 1: continue
        line = line.strip().split()
        if line[-1] == "%":
            bp = line[-4]
            cnt = line[-5]
            name = '_'.join(line[0:-5])
            if isnumber(cnt):
                if name not in dataArray:
                    keyorder.append(name)
                    dataArray[name] = [int(cnt), int(bp)]
                    continue
                dataArray[name][0] += int(cnt)
                dataArray[name][1] += int(bp)
    return dataArray, keyorder


def isnumber(v):
    try:
        a = int(v)
        return True
    except:
        try:
            a = float(v)
            return True
        except:
            return False
    

if __name__ == "__main__":
    folder = sys.argv[1]
    dirs = [f for f in os.listdir(folder) if os.path.isdir(f)]
    flist = [os.path.join(folder, fld, f) for fld in dirs for f in os.listdir(os.path.join(folder, fld)) if ".tbl" in f]
    outArray = {}
    keyorder = []
    for tblname in flist:
        outArray, keyorder = parse_tbl(tblname, outArray, keyorder)
    print("{}\t{}\t{}".format("Class", "Count", "Bp"))
    for key in keyorder:
        cnt, bp = outArray.get(key)
        print("{}\t{}\t{}".format(key, cnt, bp))
