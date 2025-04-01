import sys
import os
from os.path import exists
import filecmp
import copy
import json
import ntpath
_KB_ = 3.1668090406215046e-06
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

def getMolBuildFile(sysData, basisStr, outFname):
    with open(outFname, "w") as f:
        print('''mol = gto.Mole()''', file = f)
        print('''mol.verbose = 4''', file = f)
        print('''mol.atom="""''', file = f)
        geom = sysData["geom"]
        charge = sysData["charge"]
        spin = sysData["mult"] - 1
        for g in geom:
            line = ""
            for v in g:
                line += v + " "

            line = line.strip()
            print(line, file = f)

        print('''"""''', file = f)
        print("mol.charge = " + str(charge), file = f)
        print("mol.spin  = " + str(spin), file = f)
        print("mol.basis = " + "'" + basisStr + "'", file = f)
        print("mol.build()", file = f)
        if spin == 0:
            print("mfl = dft.RKS(mol)", file = f)
        else:
            print("mfl = dft.UKS(mol)", file = f)

        print("\n", file = f)


def getPreSCFFile(sysData, xcStr = "PBE", method = 'diis', convTol=1e-6, outFname="prescfTmp"):
    with open(outFname, "w") as f:
        spin = sysData["mult"] - 1
        if xcStr in ["PBE", "NNGGA", "SCAN", "R2SCAN"]:
            print("mfl.xc = 'PBE,PBE'", file = f)

        if xcStr in ["PW92", "NNLDA"]:
            print("mfl.xc = 'SPW92'", file = f)

        if method.lower() not in ['diis', 'newton']:
            raise ValueError(f'''Invalid method {method} passed. Valid methods: 'diis', 'newton' ''')

        if method.lower() == 'newton':
            print("mfl = mfl.newton()", file = f)

        print("mfl.max_cycle = 100", file = f)
        print(f"mfl.conv_tol = {convTol}", file = f)
        print("mfl.kernel()", file = f)
        print("dm = mfl.make_rdm1()", file = f)

        print("\n", file = f)
        if spin == 0:
            print("mfl = dft.RKS(mol)", file = f)
        else:
            print("mfl = dft.UKS(mol)", file = f)

        print("mfl.init_guess = dm", file = f)


def getSCFFile(sysData, xcStr = "PBE", convTol=1e-7, occupancy=None, ptcFile = None, tol=1e-8, sthres=0.0, printDM=True, postfix="tmp", outFname="scfTmp"):
    with open(outFname, "w") as f:
        spin = sysData["mult"] - 1
        if xcStr == "PBE":
            print("mfl.xc = 'PBE,PBE'", file = f)

        elif xcStr == "NNLDA":
            print("tol = " + str(tol), file = f)
            print('''nnlda = NNLDA("''' + ptcFile + '''", tol)''', file = f)
            print('''mfl = mfl.define_xc_(eval_xc_lda, 'LDA')''', file = f)

        elif xcStr == "NNGGA":
            print("tol = " + str(tol), file = f)
            print("sthres = " + str(sthres), file = f)
            print('''nngga = NNGGA("''' + ptcFile + '''", tol, sthres)''', file = f)
            print('''mfl = mfl.define_xc_(eval_xc_gga, 'GGA')''', file = f)

        elif xcStr == "PW92":
            print("mfl.xc = 'SPW92'", file = f)

        elif xcStr == "R2SCAN":
            print("mfl.xc = 'R2SCAN'", file = f)

        else:
            raise Exception("Invalid xcStr")

        print('''mfl.max_cycle = 100''', file = f)
        print(f'''mfl.conv_tol = {convTol}''', file = f)
        if occupancy:
            if occupancy["type"] is not None:
                otype = occupancy["type"].lower()
                if otype not in ["frac", "fermi"]:
                    raise ValueError(f"""Invalid occupancy type {otype} provide. Valid types: "frac" and "fermi" """)

                if otype == "frac":
                    print('''mfl = scf.addons.frac_occ(mfl)''', file = f)

                if otype == "fermi":
                    temp = float(occupancy["T"])
                    sigma = _KB_*temp
                    print(f'''mfl = scf.addons.smearing_(mfl, sigma={sigma}, method='fermi')''', file = f)

        print('''mfl.kernel()''', file = f)
        print('''dm = mfl.make_rdm1()''', file = f)
        if spin == 0:
            dmFname = "DM_" + postfix
            printStr = f'''"{dmFname}"'''
            print(f'''np.savetxt({printStr}, dm)''', file = f)

        else:
            dmFname = "DM0_" + postfix
            printStr = f'''"{dmFname}"'''
            print(f'''np.savetxt({printStr}, dm[0])''', file = f)
            dmFname = "DM1_" + postfix
            printStr = f'''"{dmFname}"'''
            print(f'''np.savetxt({printStr}, dm[1])''', file = f)

        print("\n", file = f)


def createScripts(systems, sysFilePostfix, header, rootDir, basisStr, xcStr, convTol = 1e-7, occupancy = None, useGuess = False, printDM = True, ptcFile = None, tol=1e-8, sthres=0.0):
    for s in systems:
        words = s.split(":")
        subset = words[0]
        sysName = words[1]
        sysDir = os.path.join(rootDir, subset, sysName)
        if not os.path.exists(sysDir):
            os.makedirs(sysDir)


        scriptFilename = sysName + sysFilePostfix + ".py"
        scriptPath = os.path.join(rootDir, subset, sysName, scriptFilename)
        scriptF = open(scriptPath, "w")
        # add the header to script
        os.system("cat " + header + " > " + scriptPath)
        molTmpFname = "molTmp"
        preSCFTmpFname = "preSCFTmp"
        scfTmpFname = "scfTmp"
        getMolBuildFile(systems[s], basisStr, molTmpFname)
        os.system("cat " + molTmpFname + " >> " + scriptPath)
        if useGuess["flag"]:
            getPreSCFFile(systems[s], useGuess["xcStr"], useGuess["method"], useGuess["convTol"], preSCFTmpFname)
            os.system("cat " + preSCFTmpFname + " >> " + scriptPath)

        getSCFFile(systems[s], xcStr, convTol, occupancy, ptcFile, tol, sthres, printDM, sysFilePostfix, scfTmpFname)
        os.system("cat " + scfTmpFname + " >> " + scriptPath)



def run(inp, stdoutF = sys.stdout):
    sysJSON = inp["sysJSON"]
    geomJSON = inp["geomJSON"]
    rootDir = os.path.abspath(inp["rootDir"])
    sysFilePostfix = inp["sysFilePostfix"]
    header = inp["header"]
    outFname = inp["outFile"]
    useGuess = inp["useGuess"]
    xc = inp["xc"]
    xcType = xc["type"]
    ptcFile = xc["ptcPath"]
    tol = xc["tol"]
    sthres = xc["sthres"]
    basisStr = inp["basis"]
    convTol = inp["convTol"]
    occupancy = inp["occupancy"]
    printDM = inp["printDM"]
    startID = inp["sysRange"]["start"]
    endID = inp["sysRange"]["end"]
    modulesAndEnv = inp["modulesAndEnv"]
    slurm = inp["slurm"]
    jobscript = os.path.abspath(inp["jobscript"]+ "_" + str(startID) + "_" + str(endID))

    print("start and end ID", startID, endID)
    # get systems based on start and end IDs
    systems = getSysData(sysJSON, geomJSON, startID, endID)
    createScripts(systems, sysFilePostfix, header, rootDir, basisStr, xcType, convTol, occupancy, useGuess, printDM, ptcFile, tol, sthres)

    print("Running for the following systems", file = stdoutF)
    print(list(systems.keys()), file = stdoutF)
    stdoutF.flush()

    nsys = len(systems)
    # if joscript exists create a new one using a number at the end
    if os.path.exists(jobscript):
        count = 2
        found = True
        while found:
            tmpname = jobscript + "_" + str(count)
            if os.path.exists(tmpname):
                count +=1
            else:
                jobscript = tmpname
                break

    with open(jobscript, "w") as f:
        print("#!/bin/sh", file = f)
        jobname = str(slurm["jobname"]) + "_" + str(startID) + "_" + str(endID)
        print("#SBATCH -A " + str(slurm["account"]), file = f)
        print("#SBATCH -J " + jobname, file = f)
        print("#SBATCH -o " + jobname+".out.%j", file = f)
        print("#SBATCH -e " + jobname+".err.%j", file = f)
        if slurm["partition"] is not None and slurm["partition"] != "":
            print("#SBATCH -p " + str(slurm["partition"]), file = f)

        if slurm["queue"] is not None and slurm["queue"] != "":
            print("#SBATCH -q " + str(slurm["queue"]), file = f)

        print("#SBATCH -t " + str(slurm["t"]), file = f)
        print("#SBATCH --nodes=" + str(slurm["nodesPerSys"]*nsys), file = f)
        for e in slurm["extras"]:
            print("#SBATCH " + e, file=f)

        print("#SBATCH --mail-type=BEGIN,END", file = f)
        print("#SBATCH --mail-user=" + str(slurm["email"]), file = f)
        print("\n", file = f)

        with open(modulesAndEnv) as f2:
            lines = f2.readlines()
            for line in lines:
                print(line.strip(), file = f)

        print("\n", file = f)
        print("export OMP_NUM_THREADS="+str(slurm["threads"]), file=f)
        for s in systems:
            words = s.split(":")
            subset = words[0]
            sysName = words[1]
            sysDir = os.path.join(rootDir, subset, sysName)
            ptcFilePath = os.path.join(rootDir, ptcFile)
            # copy ptc file to sysDir
            os.system("cp " + ptcFile + " " + sysDir)

            pyscfScript = sysName + sysFilePostfix + ".py"
            print("cd " + sysDir, file = f)
            nodesPerSys = "--nodes=" + str(slurm["nodesPerSys"])
            tasksPerNode = "--ntasks-per-node=" + str(slurm["tasksPerNode"])
            cpusPerTask = "--cpus-per-task=" + str(slurm["cpusPerTask"])
            threads = "OMP_NUM_THREADS=" + str(slurm["threads"])
            #srunLine =  threads + " srun -u " + nodesPerSys + " python " +\
            #           pyscfScript + " &> " + outFname + " &"
            srunLine = "srun -u " + nodesPerSys + " " + tasksPerNode + " " + cpusPerTask + " python " + pyscfScript + " &> " + outFname + " &"
            print(srunLine, file = f)

        print("wait", file = f)

    os.system("sbatch " + jobscript)

if __name__ == "__main__":
    correctUse = True
    usageMsg = """Usage: `python CreatePyscfScript.py inp.json [stdout_file]`,"""\
               """ where inp.json contains various input parameters,"""\
               """ and stdout_file is an optional filename to redirect the stdout."""\

    if len(sys.argv) not in [2,3,4,5]:
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

    # modify things if start and end ids are given
    if len(sys.argv) > 3:
        inp["sysRange"]["start"] = int(sys.argv[2])
        inp["sysRange"]["end"] = int(sys.argv[3])

    stdOutF = sys.stdout
    if inp["stdOut"] is not None:
        stdOutF = open(inp["stdOut"], "w")

    run(inp, stdOutF)

    stdOutF.close()
