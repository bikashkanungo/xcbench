import json

data = {}
nreactions = {}
avgE = {}
with open('GMTKN55.json') as f:
    data = json.load(f)

for subset in data:
    avgE[subset] = 0.0
    nreactions[subset] = len(data[subset])
    for index in data[subset]:
        avgE[subset] += abs(data[subset][index]["Ref."])/nreactions[subset]

nsubsets = len(data)
avgEAll = 0.0
nreactionsAll = 0
for subset in data:
    avgEAll += avgE[subset]/nsubsets
    nreactionsAll += nreactions[subset]

print('nsubsets', nsubsets)
print('nreactionsAll', nreactionsAll)
print('avgEAll', avgEAll)

subsetsToFind = ['ALK8', 'ALKBDE10', 'BHDIV10', 'HEAVYSB11', 'G2RC']
for s in subsetsToFind:
    w = avgEAll/nreactionsAll*(nreactions[s]/avgE[s])
    print(s, "Weight", w, 'avgE', avgE[s])
