import os
import json
import copy
import sys

if __name__ == "__main__":
    correctUse = True
    usageMsg = """Usage: `python GetMeanErrs.py GMTKN55.json Functionals.json` """\
               """, where GMTKN55.json is the json file created using Collate.py and manually updated for """\
               """any inconsistency/missing data (run Collate.py and follow instructions). Functionals.json """\
               """should contain all the functionals to be considered."""

    if len(sys.argv) != 3:
        if len(sys.argv) > 1:
            if sys.argv[1] == "--help" or sys.argv[1] == "-h":
                print(usageMsg)
                sys.exit(0)

        else:
            correctUse = False

    else:
        if ".json" not in sys.argv[1] or ".json" not in sys.argv[2]:
            correctUse = False

    if not correctUse:
        print("Incorrect arguments passed. " + usageMsg)
        raise RuntimeError

    dataJSON = str(sys.argv[1])
    funcJSON = str(sys.argv[2])
    outFname = "MeanErr.json"
    data = {}
    with open(dataJSON) as f:
        data = json.load(f)

    funcs = []
    with open(funcJSON) as f:
        tmp = json.load(f)
        funcs = [x for x in tmp]

    err = {}
    types = ["without", "D3(0)", "D3(BJ)"]
    for subset in data:
        err[subset] = {}
        for func in funcs:
            err[subset][func] = {}
            for t in types:
                err[subset][func][t] = []

    for subset in data:
        for systemID in data[subset]:
            sysData = data[subset][systemID]
            ref = sysData["Ref."]
            for func in funcs:
                if func in sysData:
                    for t in types:
                        if t in sysData[func] and sysData[func][t] is not None:
                            err[subset][func][t].append(sysData[func][t]-ref)
    meanErr = {}
    for subset in err:
        meanErr[subset] = {}
        for func in funcs:
            meanErr[subset][func] = {}
            for t in types:
                meanErr[subset][func][t] = {"ME": None, "MAE": None}
                N = len(err[subset][func][t])
                if N > 0:
                    me = 0.0
                    mae = 0.0
                    for v in err[subset][func][t]:
                        me += v
                        mae += abs(v)

                    me /= N
                    mae /= N
                    meanErr[subset][func][t] = {"ME": me, "MAE": mae}


    with open(outFname, "w") as f:
        json.dump(meanErr, f, indent = 2)
