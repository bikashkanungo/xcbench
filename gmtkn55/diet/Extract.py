import json
import os

inFname = "combo_150.json"
outFname = "DietGMTKN55.json"
data = {}
with open(inFname, "r") as f:
    data_ = json.load(f)
    N = 0
    NIDs = 0
    species = set()
    for k in data_.keys():
        pos1 = k.rfind('-')
        subset = k[0:pos1]
        indexStr = k[pos1+1:]
        index = int(indexStr)
        wt = data_[k]["Weight"]
        NIDs += 1
        N += len(data_[k]["Species"])
        for s in data_[k]["Species"]:
            species.add(s["ID"])

        if subset not in data.keys():
            data[subset] = {}
            data[subset]["Weight"] = wt
            data[subset]["IDs"] = [index]

        else:
            data[subset]["IDs"].append(index)

    print("NSpecies", N)
    print("NSpecies2", len(species))
    print("NIDs", NIDs)
with open(outFname, "w") as f:
   json.dump(data, f, indent=2)

