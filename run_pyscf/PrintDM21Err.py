import json

DM21ErrFname = "DM21GMTKN55.json"
weightsFname = "weights.json"

dm21Err = {}
with open(DM21ErrFname) as f:
    dm21Err = json.load(f)


weights = {}
with open(weightsFname) as f:
    weights = json.load(f)

mae = {}
wtmad = {}
for w in weights:
    wtmad[w] = 0.0
    for s in weights[w]:
        mae[s] = dm21Err[s]
        wtmad[w] += weights[w][s]*mae[s]

print(mae)
print(wtmad)

