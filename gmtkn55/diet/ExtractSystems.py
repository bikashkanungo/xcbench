import json

fullFname = "GMTKN55.json"
dietFname = "DietGMTKN55.json"
outFname = "DietGMTKN55Systems.json"

full = {}
diet = {}
with open(fullFname) as f:
    full = json.load(f)

with open(dietFname) as f:
    diet = json.load(f)

systems = {}
for subset in diet:
    sysSet = set()
    for i in diet[subset]["IDs"]:
        sysNames = full[subset][str(i)]["Systems"]
        sysSet.update(sysNames)

    systems[subset] = list(sysSet)

with open(outFname, "w") as f:
    json.dump(systems, f, indent = 2)




