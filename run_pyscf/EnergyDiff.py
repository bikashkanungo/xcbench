import sys
import json
import os
import numpy as np
from decimal import Decimal, getcontext

HARTREE_TO_KCAL_PER_MOL= 627.5096080305927
DELTA_E_TOL = 1e-4
RMSE_TOL=5e-4
NLAST = 10

def parseEnergy(sysNames, rootoutfilename, rmseTol=RMSE_TOL, NLast=NLAST):
    rootdir = os.getcwd()
    energies = {}
    for sys in sysNames:
        words = sys.split(":")
        subset = words[0]
        sysName = words[1]
        sysDir = os.path.join(rootdir, subset, sysName)
        outfilename = os.path.join(sysDir,rootoutfilename)
        #print('parsing:', outfilename)
        out = open(outfilename, "r")
        outLines = out.readlines()
        converged = False
        energy = None
        with open(outfilename) as f:
            lines = f.readlines()
            nlines = len(lines)
            for iline in range(nlines):
                # traverse in reverse
                line = lines[nlines-1-iline]
                if "converged SCF energy" in line:
                    words = line.strip().split("=")
                    energy = float(words[1].split()[0])
                    converged = True
                    break

        if not converged and rmseTol is not None:
            count = 0
            EVals = []
            with open(outfilename) as f:
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
                energy = mean
        
        if converged:
            energies[sys] = energy*HARTREE_TO_KCAL_PER_MOL
    
    return energies


def getSysNames(dietReactions, allData, subsets):
    sysNames = set()
    for subset in subsets:
        reactionIds = dietReactions[subset]["IDs"]
        subsetData = allData[subset]
        for rId in reactionIds:
            rSysNames = subsetData[str(rId)]["Systems"]
            for s in rSysNames:
                sysNames.add(subset + ":" + s)

    return list(sysNames)


def getReactionsData(dietReactions, allData, subsets):
    reactions = {}
    for subset in subsets:
        reactionIds = dietReactions[subset]["IDs"]
        subsetData = allData[subset]
        reactions[subset] = {}
        for rId in reactionIds:
            reactions[subset][rId] = subsetData[str(rId)]

    return reactions


def main(dietReactions, allData, subsets, outFileToParsei, DFAs):
    sysNames = getSysNames(dietReactions, allData, subsets)
    reactions = getReactionsData(dietReactions, allData, subsets)
    sysEnergies = parseEnergy(sysNames, outFileToParse)
    subsetsToSkip = set()
    #subsetGroups = {}
    #subsetsGroupsToConsider = ['basic', 'barrier']
    #with open("subsetsGroups.json") as f:
    #    subsetGroups = json.load(f)

    #for group in subsetGroups:
    #    if group not in subsetsGroupsToConsider:
    #        for subset in subsetGroups[group]:
    #            subsetsToSkip.add(subset)

    for sys in sysEnergies:
        words = sys.split(":")
        subset = words[0]
        sysName = words[1]
        if sysEnergies[sys] is None:
            subsetsToSkip.add(subset)
        
    err = {}
    errNN = {}
    WTMAD0 = {}
    WTMAD0["NN"] = 0.0
    WTMAD1 = {}
    WTMAD1["NN"] = 0.0
    WTMAD2 = {}
    WTMAD2["NN"] = 0.0
    for dfa in DFAs:
        WTMAD0[dfa] = 0.0
        WTMAD1[dfa] = 0.0
        WTMAD2[dfa] = 0.0

    subsetsToConsider = []
    for subset in subsets:
        if subset not in subsetsToSkip:
            subsetsToConsider.append(subset)
  
    subsetsToRemove = ['YBDE18', 'AL2X6', 'NBPRC', 'DC13', 'ALK8']
    for subset in subsetsToRemove:
        subsetsToConsider.remove(subset)

    print("subsets to consider:", subsetsToConsider)
    nsysSelected = 0
    for subset in subsetsToConsider:
        subsetData = allData[subset]
        nsysSelected += len(subsetData)

    nsubsetsConsidered = len(subsetsToConsider)
    weights1 = {} 
    weights2 = {} 
    avgRef = {} 
    nsubsetReactions = {} 
    nallReactions = 0
    allAvgRef = 0.0
    for subset in subsetsToConsider:
        err[subset] = {"NN": {"ME": 0.0, "MAE": 0.0}}
        errNN[subset] = {}
        for DFA in DFAs:
            err[subset][DFA] = {"ME": 0.0, "MAE": 0.0}

        reactionsSubset = reactions[subset]
        nreactions = len(reactionsSubset)
        nsubsetReactions[subset] = nreactions
        nallReactions += nreactions
        sumRef = 0.0
        for rId in reactionsSubset:
            rSysNames = reactionsSubset[rId]["Systems"]
            rStoch = reactionsSubset[rId]["Stoichiometry"]
            ref = reactionsSubset[rId]["Ref."]
            sumRef += abs(ref)
            nnVal = 0.0
            for i in range(len(rSysNames)):
                sysName = subset + ":" + rSysNames[i]
                energy = sysEnergies[sysName]
                nnVal += energy*rStoch[i]

            errNN[subset][rId] = nnVal - ref
            err[subset]["NN"]["ME"] += (nnVal - ref)/nreactions
            err[subset]["NN"]["MAE"] += abs(nnVal - ref)/nreactions
            for DFA in DFAs:
                dfaVal = reactionsSubset[rId][DFA]["without"]
                err[subset][DFA]["ME"] += (dfaVal - ref)/nreactions
                err[subset][DFA]["MAE"] += abs(dfaVal - ref)/nreactions

        avgRef[subset] = sumRef/nreactions
        allAvgRef += avgRef[subset]/nsubsetsConsidered

    
    for subset in subsetsToConsider:
        if avgRef[subset] > 75.0:
            weights1[subset] = 0.1/nsubsetsConsidered
        
        elif avgRef[subset] < 7.5:
            weights1[subset] = 10.0/nsubsetsConsidered

        else:
            weights1[subset] = 1.0/nsubsetsConsidered
    
        weights2[subset] = (allAvgRef/nallReactions)*(nsubsetReactions[subset]/avgRef[subset])

    with open("weights.json", "w") as fw:
        tmp = {"w1": weights1, "w2": weights2}
        json.dump(tmp, fw, indent=2)

    for s in err:
        print(s, end =  " ")
        #w0 = dietReactions[s]["Weight"]
        w1 = weights1[s]
        w2 = weights2[s]
        for dfa in err[s]:
            #WTMAD0[dfa] += w0*err[s][dfa]["MAE"]/nsubsetsConsidered
            WTMAD1[dfa] += w1*err[s][dfa]["MAE"]
            WTMAD2[dfa] += w2*err[s][dfa]["MAE"]
            print(dfa+":" , "{0:.3f}".format(err[s][dfa]["MAE"]), end = "  ")
        
        print("Weight2:", w2, 'Nreactions:', nsubsetReactions[s], "AvgRef.:", avgRef[s])
    
    print("\n\n")
  
    
    DFAsPlusNN = DFAs + ["NN"]
    # +2 for WTMAD-2 and WTMAD-2
    errMat = np.zeros((len(DFAsPlusNN), nsubsetsConsidered+2)) 
    for idfa, dfa in enumerate(DFAsPlusNN):
        for isub, s in enumerate(err):
            errMat[idfa, isub] = err[s][dfa]["MAE"]


    for idfa, dfa in enumerate(DFAsPlusNN):
        errMat[idfa, nsubsetsConsidered] = WTMAD1[dfa]
        errMat[idfa, nsubsetsConsidered+1] = WTMAD2[dfa]
    
    with open("errMat", "w") as ferrMat:
        cols = list(err.keys()) + ["WTMAD-1", "WTMAD-2"]
        print("Cols", cols, file = ferrMat)
        print("Rows", DFAsPlusNN, file= ferrMat)
        print(errMat.tolist(), file=ferrMat)
    
    print("WTMAD0", WTMAD0)
    print("WTMAD1", WTMAD1)
    print("WTMAD2", WTMAD2)
    
    print("\nNN worse than PBE")
    print("-----------------")
    for s in err:
        if err[s]["NN"]["MAE"] > err[s]["PBE"]["MAE"]:
            print(s, "NN:", err[s]["NN"]["MAE"], "PBE:", err[s]["PBE"]["MAE"])

    print("\nNN better than PBE")
    print("-----------------")
    for s in err:
        if err[s]["NN"]["MAE"] < err[s]["PBE"]["MAE"]:
            print(s, "NN:", err[s]["NN"]["MAE"], "PBE:", err[s]["PBE"]["MAE"])
    
    print("\nNN worse than SCAN")
    print("-----------------")
    for s in err:
        if err[s]["NN"]["MAE"] > err[s]["SCAN"]["MAE"]:
            print(s, "NN:", err[s]["NN"]["MAE"], "SCAN:", err[s]["SCAN"]["MAE"])
    
    print("\nNN better than SCAN")
    print("-----------------")
    for s in err:
        if err[s]["NN"]["MAE"] < err[s]["SCAN"]["MAE"]:
            print(s, "NN:", err[s]["NN"]["MAE"], "SCAN:", err[s]["SCAN"]["MAE"])

    #print("\n\n")
    #for s in errNN:
    #    print(s,errNN[s])

    

if __name__ == "__main__":
    dietJSON = str(sys.argv[1])
    allDataJSON = str(sys.argv[2])
    subsetsFile = str(sys.argv[3])
    outFileToParse = str(sys.argv[4])
    DFAs = ["PBE", "SCAN", "B3LYP", "wB97X-V"]

    dietReactions = None
    with open(dietJSON) as f:
        dietReactions = json.load(f)

    allData = None
    with open(allDataJSON) as f:
        allData = json.load(f)

    subsets = []
    with open(subsetsFile) as f:
        lines = f.readlines()
        for line in lines:
            if line:
                words = line.strip().split()
                subsets.append(words[0])

    main(dietReactions, allData, subsets, outFileToParse, DFAs)
