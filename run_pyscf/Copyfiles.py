import json
import os

rootdir = "./"
inpJSON = "UnconvergedSys_out_PW92_def2_qzvp.json"
srcFname = "out"
dstFname =  
inp = None
with open(inpJSON) as f:
    inp = json.load(f)

for subset in inp:
    for sys in inp[subset]:
        path = os.path.join(rootdir, subset, sys)
        os.chdir(path)
        os.system("cp " + srcFname + " " + dstFname)

