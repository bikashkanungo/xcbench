import json
import sys


if __name__ == "__main__":
    inp = {}
    inpFname = str(sys.argv[1])
    with open(inpFname) as f:
        inp = json.load(f)


    allJSONFname = inp["allJSON"]
    outSysJSONFname = inp["outSysJSON"]
    subsetsGroupsJSONFname = inp["subsetsGroupsJSON"]
    groups = inp["groups"]
    subsetsGroups = {}
    allData = {}
    with open(allJSONFname) as f:
        allData = json.load(f)

    with open(subsetsGroupsJSONFname) as f:
        subsetsGroups = json.load(f)

    systems = {}
    for group in groups:
        subsets = subsetsGroups[group]
        for subset in subsets:
            systems[subset] = set()
            subsetData = allData[subset]
            for index in subsetData:
                for s in subsetData[index]["Systems"]:
                    systems[subset].add(s)
    
    outDict = {}
    for subset in systems:
        outDict[subset] = list(systems[subset])
    
    with open(outSysJSONFname, "w") as f:
        json.dump(outDict, f, indent = 2)

  
