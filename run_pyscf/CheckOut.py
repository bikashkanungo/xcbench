import sys
import os
from os.path import exists
import filecmp
import copy
import json
import ntpath
import numpy as np
    
def getUniqueFname(prefix, rootDir):
    # if prefix exists, created a new one with name prefix_<N>,
    #where <N> is the lowest positive integer such that jobscript_<N> does not exist
    fpath = os.path.join(rootDir, prefix)
    if os.path.exists(fpath):
        created = False
        count = 1
        while not created:
            prefixNew = prefix +"_" + str(count)
            fpath = os.path.join(rootDir, prefixNew)
            if not os.path.exists(fpath):
                created = True

            else:
                count += 1

        return prefixNew

    else:
        return prefix

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
    usageMsg = """Usage: `python Checkout.py inp.json [RMSE_tol] [N_last]`,"""\
               """ where inp.json contains various input parameters. """\
               """This script checks if the string "converged SCF energy" is found in the """\
               """pyscf output file. In not, it can check the energy values from the """\
               """last N (given by N_last) SCF iterations and deems it converge if the RMSE of """\
               """those N_last values is less than the RMSE_tol."""

    if len(sys.argv) not in [2,4]:
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
    rmseTol = None 
    NLast = None 
    if len(sys.argv) == 4:
        rmseTol = float(sys.argv[2])
        NLast = int(sys.argv[3])
        print(f"Using RMSE_tol = {rmseTol} and N_last = {NLast} ")

    inp = {}
    with open(inpJSON) as f:
        inp = json.load(f)
    
    sysJSON = inp["sysJSON"]
    geomJSON = inp["geomJSON"]
    startID = inp["sysRange"]["start"]
    endID = inp["sysRange"]["end"]
    outFname = inp["outFile"]
    rootDir = os.path.abspath(inp["rootDir"])
    # get systems based on start and end IDs
    systems = getSysData(sysJSON, geomJSON, startID, endID)

    unconvergedSys = []
    for s in systems:
        words = s.split(":")
        subset = words[0]
        sysName = words[1]
        sysDir = os.path.join(rootDir, subset, sysName)
        outFilepath = os.path.join(sysDir, outFname)
        converged = False
        with open(outFilepath) as f:
            lines = f.readlines()
            for line in reversed(lines):
                if "converged SCF energy" in line:
                    converged = True
                    break
        
        if not converged and rmseTol is not None:
            count = 0
            EVals = []
            with open(outFilepath) as f:
                lines = f.readlines()
                for line in reversed(lines):
                    if "cycle=" in line:
                        index = line.index('E=') 
                        E = float(line[index:].split()[1])
                        EVals.append(E)
                        count += 1

                    if count >= NLast:
                        break

            EVals = np.array(EVals)
            mean = EVals.sum()/count
            rmse = ((((EVals - mean)**2.0).sum())/count)**0.5
            if rmse < rmseTol:
                converged = True
            
        if not converged:
            unconvergedSys.append(s)

    print("Unconverged systems:", unconvergedSys)

    if len(unconvergedSys) > 0:
        unconvergedDict = {}
        for s in unconvergedSys:
            words = s.split(":")
            subset = words[0]
            sysName = words[1]
            if subset not in unconvergedDict:
                unconvergedDict[subset] = [sysName]

            else:
                unconvergedDict[subset].append(sysName)

        unconvFPrefix = "UnconvergedSys_" + outFname
        unconvFPrefix = getUniqueFname(prefix = unconvFPrefix, rootDir = "./")
        unconvJSON = unconvFPrefix + ".json"
        with open(unconvJSON, "w") as ff:
            json.dump(unconvergedDict, ff, indent=2)

