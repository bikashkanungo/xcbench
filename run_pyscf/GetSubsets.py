import sys
import os
from os.path import exists
import filecmp
import copy
import json
import ntpath

def getSysData(sysJSON, geomJSON, start=None, end=None, includeChargedSys = True, chargeTol = 1e-6):
    data = None
    sysNames = []
    systems = {}
    with open(sysJSON) as f:
        tmp = json.load(f)
        for subset in tmp:
            for s in tmp[subset]:
                sysNames.append(subset + ":" + s)

    if start is None:
        start = 0

    if end is None:
        end = len(sysNames)

    # consider only the systems between start and end
    sysNames = sysNames[start:end]

    # load the data for all systems
    with open(geomJSON) as f:
        data = json.load(f)

    for s in sysNames:
        words = s.split(":")
        subset = words[0]
        sysName = words[1]
        systems[s] = data[subset][sysName]

    return systems

if __name__ == "__main__":
    correctUse = True
    usageMsg = """Usage: `python GetSys.py inp.json [stdout_file]`,"""\
               """ where inp.json contains various input parameters,"""\
               """ and stdout_file is an optional filename to redirect the stdout."""\

    if len(sys.argv) not in [2,3]:
        correctUse = False

    else:
        if sys.argv[1] == "--help" or sys.argv[1] == "-h":
            print(usageMsg)
            sys.exit(0)

        else:
            if ".json" not in sys.argv[1]:
                correctUse  = False

    if not correctUse:
        print("Incorrect arguments. """ + usageMsg)
        raise RuntimeError

    inpJSON = str(sys.argv[1])
    if len(sys.argv) > 2:
        stdoutF = open(sys.argv[2], "w")

    else:
        stdoutF = sys.stdout

    inp = {}
    with open(inpJSON) as f:
        inp = json.load(f)
    
    sysJSON = inp["sysJSON"]
    geomJSON = inp["geomJSON"]
    startID = inp["sysRange"]["start"]
    endID = inp["sysRange"]["end"]

    subsets = set()
    # get systems based on start and end IDs
    systems = getSysData(sysJSON, geomJSON, startID, endID)

    for s in systems:
        words = s.split(":")
        subsets.add(words[0])

    print("Subsets between system range[", startID, ",", endID, ") :")
    for subset in subsets:
        print(subset)
