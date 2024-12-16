import sys
import os
from os.path import exists
import filecmp
import copy
import json
import ntpath
import PeriodicTable as ptab

__ANGS_TO_BOHR__ = 1.8897259885789

def getSysData(sysJSON, geomJSON, start=None, end=None):
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

def createCoords(systems, rootDir, coordsFname, domainFname, vaccuum=25.0, base = 5):
    for s in systems:
        words = s.split(":")
        subset = words[0]
        sysName = words[1]
        sysData = systems[s]
        mag = 0.0
        if sysData["mult"] > 1:
            mag = -0.05

        path = os.path.join(rootDir, subset, sysName)
        if not os.path.exists(path):
            os.makedirs(path)

        os.chdir(path)
        minCoord = [1e12 for i in range(3)]
        maxCoord = [-1e12 for i in range(3)]
        with open(coordsFname, "w") as f:
            geom = sysData["geom"]
            for atom in geom:
                sym = atom[0]
                # safeguard for proper symbols.
                # Example: convert "CL" to "Cl" for chlorine
                # or convert "c" to "C" for carbon
                if len(sym) == 1:
                    sym = sym.upper()
                if len(sym) > 1:
                    sym = (sym[0]).upper() + (sym[1:]).lower()

                Z = ptab.getZ(sym)
                coords = [float(x)*__ANGS_TO_BOHR__ for x in atom[1:]]
                print(Z, Z, coords[0], coords[1], coords[2], mag, file = f)
                for i in range(3):
                    minCoord[i] = min(minCoord[i], coords[i])
                    maxCoord[i] = max(maxCoord[i], coords[i])

        with open(domainFname, "w") as f:
            for i in range(3):
                dv = [0.0, 0.0, 0.0]
                dv[i] = 2*(maxCoord[i] - minCoord[i] + vaccuum)
                # roundoff to the base
                dv[i] = base * round(dv[i]/base)
                print(dv[0], dv[1], dv[2], file = f)

    # return back to root dir
    os.chdir(rootDir)

def getMeshParams(Z):
    meshParams = {'feorder': 4, 'meshAroundAtom' : 0.4,  'atomBallRadius': 1.5, 'meshAtAtom': 0.05}
    meshParams['feorderElectrostatic'] = meshParams['feorder'] + 2
    if Z > 8 and Z < 13:
        meshParams['feorder'] = 5
        meshParams['feorderElectrostatic'] = meshParams['feorder'] + 2
        meshParams['meshAroundAtom'] = 0.4
        meshParams['atomBallRadius'] = 1.5
        meshParams['meshAtAtom'] = 0.02

    if Z >= 13:
        meshParams['feorder'] = 6
        meshParams['feorderElectrostatic'] = meshParams['feorder'] + 2
        meshParams['meshAroundAtom'] = 0.4
        meshParams['atomBallRadius'] = 1.5
        meshParams['meshAtAtom'] = 0.02

    return meshParams

def getAtomicNumbers(sysData):
    geom = sysData["geom"]
    atomicNums = []
    for atom in geom:
        sym = atom[0]
        # safeguard for proper symbols.
        # Example: convert "CL" to "Cl" for chlorine
        # or convert "c" to "C" for carbon
        if len(sym) == 1:
            sym = sym.upper()

        if len(sym) > 1:
            sym = (sym[0]).upper() + (sym[1:]).lower()

        Z = ptab.getZ(sym)
        atomicNums.append(Z)

    return atomicNums

def writeParamsFile(keysToWordsAndVals, paramsLines, filePath):
    keys  = list(keysToWordsAndVals.keys())
    nkeys = len(keys)
    keywords = []
    vals = []
    for k in keys:
        keywords.append(keysToWordsAndVals[k][0])
        vals.append(keysToWordsAndVals[k][1])

    with open(filePath, "w") as f:
        for paramLine in paramsLines:
            lineToPrint = paramLine.rstrip()
            equalToPos = paramLine.find("=")
            if equalToPos != -1:
                keyword = (paramLine[:equalToPos]).strip()
                for i in range(nkeys):
                    if keyword == keywords[i]:
                        lineToPrint = paramLine[:equalToPos] + " = " + str(vals[i])

            print(lineToPrint, file = f)


def createModelFile(xcType, modelPath, ptcPath, sysDir):
    modelFname = ntpath.basename(modelPath)
    ptcFname = ntpath.basename(ptcPath)
    os.system("cp " + ptcPath + " " + sysDir)
    inpF = open(modelPath, 'r')
    outF = open(os.path.join(sysDir, modelFname), 'w')
    lines = inpF.readlines()
    for line in lines:
        line_ = line
        words = line.strip().split("=")
        if words[0] == "PTC_FILE":
            ptcFound = words[1]
            if ptcFound != ptcPath and pctFound != ptcFname:
                line_ = words[0] + " = " + ptcFname

        print(line_, file = outF)

    inpF.close()
    outF.close()

def createFiles(inp):
    inpParamsFname = inp["inpParamsFile"]
    outParamsFname = inp["outParamsFile"]
    sysJSON = inp["sysJSON"]
    geomJSON = inp["geomJSON"]
    rootDir = os.path.abspath(inp["rootDir"])
    coordsFname = inp["outCoordsFile"]
    domainFname = inp["outDomainFile"]
    vaccuumBuffer = inp["vaccuumBuffer"]
    xc = inp["xc"]
    xcType = xc["type"]
    modelPath = xc["modelPath"]
    ptcPath = xc["ptcPath"]
    startID = inp["sysRange"]["start"]
    endID = inp["sysRange"]["end"]
    dispersionType = inp["dispersion"]
    d3DampingType = inp["D3Damping"]

    paramsLines = None
    with open(inpParamsFname, 'r') as f:
        paramsLines = f.readlines()


    paramsNewFilename = "parameters.prm"

    keysToWordsAndDefaultVals = {"natoms" : ["set NATOMS",0],
                                 "ntypes" : ["set NATOM TYPES", 0],
                                 "coordsFile": ["set ATOMIC COORDINATES FILE", "coordinates.inp"],
                                 "domainFile": ["set DOMAIN VECTORS FILE", "domainVectors.inp"],
                                 "maxIter" : ["set MAXIMUM ITERATIONS", 80],
                                 "temp" : ["set TEMPERATURE", 100],
                                 "scfTol" : ["set TOLERANCE", 1e-4],
                                 "chebTol": ["set CHEBYSHEV FILTER TOLERANCE", 1e-2],
                                 "meshAroundAtom" : ["set MESH SIZE AROUND ATOM", 0.2],
                                 "atomBallRadius" : ["set ATOM BALL RADIUS", 2.5],
                                 "meshAtAtom" : ["set MESH SIZE AT ATOM",0.02],
                                 "feorder" : ["set POLYNOMIAL ORDER", 4],
                                 "feorderElectrostatic": ["set POLYNOMIAL ORDER ELECTROSTATICS", 6],
                                 "spin" : ["set SPIN POLARIZATION", 0],
                                 "charge": ["set NET CHARGE", 0.0],
                                 "dispersion": ["set DISPERSION CORRECTION TYPE", 0],
                                 "D3DampingType": ["set D3 DAMPING TYPE", "0"],
                                 "xc": ["set EXCHANGE CORRELATION TYPE", "GGA-PBE"],
                                 "xcModelFile": ["set MODEL XC INPUT FILE", " "]}

    keys = list(keysToWordsAndDefaultVals)
    keywords = [v[0] for v in keysToWordsAndDefaultVals.values()]
    nkeys = len(keys)

    # sanity check for existence of files
    if 'ML' in xcType:
        if not os.path.exists(modelPath):
            raise Exception("Model file " + modelPath + " not found")

        if not os.path.exists(ptcPath):
            raise Exception("PTC file " + ptcPath + " not found")

    # get systems based on start and end IDs
    systems = getSysData(sysJSON, geomJSON, startID, endID)

    # create system folders (if they do not exists) and create
    # coordinates and domainVectors files in them
    createCoords(systems, rootDir, coordsFname, domainFname, vaccuumBuffer)
    for s in systems:
        words = s.split(":")
        subset = words[0]
        sysName = words[1]
        sysDir = os.path.join(rootDir, subset, sysName)

        sysData = systems[s]
        charge = sysData["charge"]
        mult = sysData["mult"]

        # copy the default key value pairs
        keysToWordsAndDefaultValsUpdated = copy.deepcopy(keysToWordsAndDefaultVals)

        # modify number of atoms and number of atom types
        atomicNums = getAtomicNumbers(sysData)
        atomicNumsSet = set(atomicNums)
        keysToWordsAndDefaultValsUpdated['natoms'][1] = len(atomicNums)
        keysToWordsAndDefaultValsUpdated['ntypes'][1] = len(atomicNumsSet)

        # modify coordinates and domain vectors file names
        keysToWordsAndDefaultValsUpdated["coordsFile"][1] = coordsFname
        keysToWordsAndDefaultValsUpdated["domainFile"][1] = domainFname

        # modify for spin polarized systems
        if mult > 1:
            keysToWordsAndDefaultValsUpdated['spin'][1] = 1
            keysToWordsAndDefaultValsUpdated['temp'][1] = 1

        # add any charge
        if abs(charge) > 1e-16:
            # In DFT-FE we assign positive charge when electrons are added
            # However, in usual databases addition of electrons is
            # assigned negative charge. Hence we take the negative of the
            # input charge
            keysToWordsAndDefaultValsUpdated['charge'][1] = -charge

        # get mesh params as per maximum atomic number.
        # modify getMeshParams() function accordingly
        maxAtomicNum = max(atomicNums)
        meshParams = getMeshParams(maxAtomicNum)
        for k in meshParams:
            keysToWordsAndDefaultValsUpdated[k][1] = meshParams[k]

        # modify xc type
        keysToWordsAndDefaultValsUpdated['xc'][1] = xcType

        # if using MLXC, copy the model.inp (modelPath) and ptc file (ptcPath)
        # to the system directory
        if 'ML' in xcType:
            createModelFile(xcType, modelPath, ptcPath, sysDir)
            modelFname = ntpath.basename(modelPath)
            keysToWordsAndDefaultValsUpdated['xcModelFile'][1] = modelFname

        # modify any dispersion correction parameters
        if dispersionType is not None:
            if dispersionType.lower() == "d3":
                keysToWordsAndDefaultValsUpdated['dispersion'][1] = 1
                if int(d3Damping) == 0:
                    keysToWordsAndDefaultValsUpdated['D3DampingType'][1] = 0

                elif d3Damping.lower() == 'bj':
                    keysToWordsAndDefaultValsUpdated['D3DampingType'][1] = 1

                else:
                    raise Exception("Invalid D3 damping type provided. Valid types: 0, BJ")

            if dispersionType.lower() == "d4":
                keysToWordsAndDefaultValsUpdated['dispersion'][1] = 2

        # write the new parameters file in the system directory
        filePath = os.path.join(sysDir, outParamsFname)
        writeParamsFile(keysToWordsAndDefaultValsUpdated, paramsLines, filePath)

def run(inp):
    outParamsFname = inp["outParamsFile"]
    sysJSON = inp["sysJSON"]
    geomJSON = inp["geomJSON"]
    rootDir = os.path.abspath(inp["rootDir"])
    coordsFname = inp["outCoordsFile"]
    domainFname = inp["outDomainFile"]
    outFname = inp["outFile"]
    xc = inp["xc"]
    xcType = xc["type"]
    modelPath = xc["modelPath"]
    ptcPath = xc["ptcPath"]
    startID = inp["sysRange"]["start"]
    endID = inp["sysRange"]["end"]
    dispersionType = inp["dispersion"]
    d3DampingType = inp["D3Damping"]
    dftfeExec = inp["dftfeExec"]
    srunParams = inp["srunParams"]

    if not os.path.exists(dftfeExec):
        raise Exception("DFT-FE exectuable " + dftfeExec + " not found.")

    if srunParams is None:
        srunParams = ""

    # get systems based on start and end IDs
    systems = getSysData(sysJSON, geomJSON, startID, endID)
    for s in systems:
        words = s.split(":")
        subset = words[0]
        sysName = words[1]
        sysDir = os.path.join(rootDir, subset, sysName)

        os.chdir(sysDir)

        modelFname = ntpath.basename(modelPath)
        ptcFname = ntpath.basename(ptcPath)
        # sanity checks for files
        files = [outParamsFname, coordsFname, domainFname, modelFname, ptcFname]
        for name in files:
            if not os.path.exists(name):
                msg = "In subset " + subset + " and system " + sysName + ", file " + name + " not found."
                raise Exception(msg)

        os.system("srun " + srunParams + " " + dftfeExec + " " + outParamsFname + " &> " + outFname)

if __name__ == "__main__":
    correctUse = True
    usageMsg = """Usage: `python CreateAndRun inp.json`, where inp.json contains various input parameters."""\
               """ Refer to README.md for details of the content of inp.json."""
    if len(sys.argv) != 2:
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
    inp = {}
    with open(inpJSON) as f:
        inp = json.load(f)

    mode = inp["mode"]
    if mode == "create":
        createFiles(inp)

    if mode == "run":
        run(inp)

    if mode == "both":
        createFiles(inp)
        run(inp)

