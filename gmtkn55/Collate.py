import os
import json
import copy
from bs4 import BeautifulSoup
import re
import sys

def isInt(s):
    try:
        int(s)
    except ValueError:
        return False
    else:
        return True

def isFloat(s):
    try:
        float(s)
    except ValueError:
        return False
    else:
        return True


def hasSpecialChar(s):
    return any(not c.isalnum() for c in s)

def getNames(jsonFname):
    d = {}
    with open(jsonFname) as f:
        d = json.load(f)

    return list(d.keys())

def processFuncHTML(fname, refData, subset, func, warningF):
    funcData = {}
    refDataTmp = {}
    valKeysFunc = ["without", "D3(0)", "D3(BJ)"]
    valKeysRef = ["Ref."]
    with open(fname) as f:
        # convert html to text and retain the newlines
        txt = BeautifulSoup(f, features="lxml").get_text("\n")

        # Consider the part after "Systems" keyword
        txt = txt[txt.find("Systems"):]
        header = txt[0:txt.find("1")]
        header = header.split()

        # The first two columns are Systems and Stochiometry
        # and the next are the values
        valKeys = header[2:]
        nvals = len(valKeys)
        valKeysColIds= {}
        for valKey in valKeys:
            valKeysColIds[valKey] = header.index(valKey) - 2 #len(header)

        # split based on double new line
        txt = txt.split("\n\n")
        # consider lines which start with an integer
        sysLines = []
        for line in txt:
            words = line.split()
            if words:
                if isInt(words[0]):
                    sysLines.append(line)

        for line in sysLines:
            words = line.split()
            # Sometimes there are stray tilde at the end,
            # possibly because of incorrect html to text conversion.
            # remove them explicitly
            if words[-1] == '~':
                words = words[:-1]

            # first word is the index of the reaction
            index = int(words[0])
            refDataSys = refData[index]
            sysNamesRef = refDataSys["Systems"]
            nsys = len(sysNamesRef)

            #consider the rest of the line
            words = words[1:]
            nwords = len(words)

            sysNames = words[0:nsys]
            stoch = []
            for w in words[nsys:2*nsys]:
                if isInt(w):
                    stoch.append(int(w))
                else:
                    stoch.append(None)

            vals = {}
            for valKey in valKeysFunc:
                vals[valKey] = None

            vals["Ref."] = refData[index]["Ref."]

            posTab = [m.start() for m in re.finditer('\t', line)]
            valWords = words[-nvals:]
            nvalsFound = len(posTab) + 1
            if nvals != nvalsFound:
                print("Inconsistent:", subset, func, index, "nvals:", nvals, "nvalsFound:", nvalsFound, file = warningF)
                founAllVals = False

            else:
                ii = 0
                for valKey in valKeys:
                    if valKey != "Ref.":
                        v = (line[posTab[ii]+1:]).split()[0]
                        vals[valKey] = float(v)
                        ii += 1

                    if valKey not in valKeysFunc and valKey != "Ref.":
                        vals["without"] = vals[valKey]
                        del vals[valKey]

                # GMTKN55 database has the difference (DFA_result-Reference).
                # We store the DFA_result by adding the Reference value to the difference
                for valKey in valKeysFunc:
                    if vals[valKey] is not None:
                        vals[valKey] += refData[index]["Ref."]

            funcData[index] = {}
            refDataTmp[index] = {}
            for valKey in valKeysFunc:
                funcData[index].update({valKey:vals[valKey]})

            refDataTmp[index] = {"Systems": sysNames, "Stoichiometry": stoch}
            for valKey in valKeysRef:
                refDataTmp[index].update({valKey:vals[valKey]})

    return funcData, refDataTmp


def processRefHTML(fname):
    refData = {}
    with open(fname) as f:
        # convert html to text and retain the newlines
        txt = BeautifulSoup(f, features="lxml").get_text("\n")

        # Consider the part after "Systems" keyword
        txt = txt[txt.find("Systems"):]
        header = txt[0:txt.find("1")]
        header = header.split()
        #The first two columns are system ID and Stochiometry
        nvals = len(header) - 2

        # split based on double new line
        txt = txt.split("\n\n")

        # consider lines which start with an integer
        sysLines = []
        for line in txt:
            words = line.split()
            if words:
                if isInt(words[0]):
                    sysLines.append(line)

        for line in sysLines:
            words = line.split()
            # first word is the index of the reaction
            index = int(words[0])

            #consider the rest of the line
            words = words[1:]
            nwords = len(words)
            nsys = (nwords - nvals) // 2
            if 2*nsys + nvals != nwords:
                raise Exception("Inconsistency in number of columns and number of systems. Line:" + line)

            sysNames = words[0:nsys]
            stoch = [int(x) for x in words[nsys:2*nsys]]
            ref = float(words[2*nsys])
            refData[index] = {"Systems": sysNames, "Stoichiometry": stoch, "Ref.": ref}

    return refData

def processSystems(geomDir, geomFname, chargeMultPath):
    systems = {}
    with open(chargeMultPath) as f:
        lines = f.readlines()
        for line in lines:
            words = line.split()
            systems[words[0]] = {'charge': int(words[1]), 'mult': int(words[2])}

    for name in systems:
        path = os.path.join(geomDir, name, geomFname)
        geom = []
        with open(path) as f:
            lines = f.readlines()
            # skip first two lines
            for line in lines[2:]:
                words = line.split()
                geom.append(words)

        systems[name].update({"geom":geom})

    return systems


def getMismatchedRefData(refData, refDataTmp, warningF):
    ids = []
    for index in refData:
        matched = True
        if index in refDataTmp:
            if refDataTmp[index] != refData[index]:
                matched  = False

        else:
            matched = False

        if not matched:
            ids.append(index)

    return ids

def main(toFetch = True):
    subsetFname = "Subsets.json"
    funcsFname = "Functionals.json"
    url = {"geom": {"prefix": "http://www.thch.uni-bonn.de/tc.old/downloads/GMTKN/GMTKN55/", "postfix":".tar"},
           "refData": {"prefix": "http://www.thch.uni-bonn.de/tc.old/downloads/GMTKN/GMTKN55/", "postfix":"ref.html"},
           "chargeMult": {"prefix": "http://www.thch.uni-bonn.de/tc.old/downloads/GMTKN/GMTKN55/CHARGE_MULTIPLICITY_", "postfix":".txt"},
           "funcData": {"prefix": "http://www.thch.uni-bonn.de/tc.old/downloads/GMTKN/GMTKN55/results/", "postfix":"/result.html"}
           }

    geomDirname = "geom"
    geomFname = "struc.xyz"
    chargeMultFname = "ChargeMult.txt"
    refHtmlFname = "ref.html"
    tmpTar = "tmp.tar"
    outGeom = "Geom.json"
    outData = "Data.json"
    warningF = open("Warning", "w")
    inconsistencyF = open("Inconsistency", "w")

    subsets = getNames(subsetFname)
    functionals = getNames(funcsFname)

    rootDir = os.getcwd()

    geomAllSubsets = {}
    dataAllSubsets = {}
    for subset in subsets:
        path = os.path.join(rootDir, subset)
        if not os.path.exists(path):
            os.mkdir(path)

        os.chdir(path)


        # get the geometry tar file and extract it
        geomDirPath = os.path.join(os.getcwd(), geomDirname)

        if toFetch:
            if not os.path.exists(geomDirPath):
                os.mkdir(geomDirPath)

            # First extract the tar files to tmp.tar
            tarUrl = url["geom"]["prefix"] + subset + url["geom"]["postfix"]
            cmd1 = "wget -O " + tmpTar + " " + tarUrl
            os.system(cmd1)

            # Next, extract the tmp.tar to the geomDirPath folder
            # We do this because the underlying name of the folder upon untaring a .tar file
            # is not known to us.
            cmd2 = "tar -xf " + tmpTar + " -C " + geomDirPath + " --strip-components=1"
            os.system(cmd2)

            # get charge and spin multiplicity data
            os.system("wget -O " + chargeMultFname + " " + url["chargeMult"]["prefix"] + subset + url["chargeMult"]["postfix"])

        #process systems in the subset
        chargeMultPath = os.path.join(os.getcwd(), chargeMultFname)
        systems = processSystems(geomDirPath, geomFname, chargeMultPath)
        geomAllSubsets[subset] = systems

        if toFetch:
            # store it as a tmp.html
            os.system("wget -O " + refHtmlFname + " " + url["refData"]["prefix"] + subset + url["refData"]["postfix"])

        htmlPath = os.path.join(os.getcwd(), refHtmlFname)
        # process reference data
        refData = processRefHTML(htmlPath)

        # store reference data
        data = {}
        for index in refData:
            data[index] = copy.deepcopy(refData[index])

        # get functional data and process it
        for ifunc, func in enumerate(functionals):
            # store it as a tmp.html
            funcHtmlFname = func + ".html"
            if toFetch:
                os.system("wget -O " + funcHtmlFname + " " + url["funcData"]["prefix"] + subset + "/" + func + url["funcData"]["postfix"])

            htmlPath = os.path.join(os.getcwd(), funcHtmlFname)
            funcData, refDataTmp = processFuncHTML(htmlPath, refData, subset, func, inconsistencyF)
            if refData != refDataTmp:
                misIds = getMismatchedRefData(refData, refDataTmp, warningF)
                print("For subset", subset + ", mismatch in reference data for functional", func,\
                        "for the follwing indices:", misIds, file = warningF)

            for index in refData:
                data[index].update({func: funcData[index]})

        dataAllSubsets[subset] = data


    # move to rootDir
    os.chdir(rootDir)
    with open(outGeom, "w") as f:
        json.dump(geomAllSubsets, f, indent = 2)

    with open(outData, "w") as f:
        json.dump(dataAllSubsets, f, indent = 2)

    warningF.close()
    inconsistencyF.close()

if __name__ == "__main__":
    toFetch = True
    correctUse = True
    usageMsg = "Usage: `python Collate.py <fetch-bool>`, """\
                """where the <fetch-bool> can be true or false (case insensitive). """\
                """If true, it fetches all the data from the GMTKN55.json database """\
                """into individual subset directories and then collated them. """\
                """If false, it assumes that the data files and been already fetched """\
                """and the collation process in performed."""

    if len(sys.argv) != 2:
        correctUse = False

    else:
        if str(sys.argv[1]) == "--help" or str(sys.argv[1]) == "-h":
            print(usageMsg)
            sys.exit(0)

        else:
            if str(sys.argv[1]).lower()[0] not in ["t", "f"]:
                correctUse = False

    if not correctUse:
        print ("""Incorrect arguments passed. """ + usageMsg)
        raise RuntimeError

    else:
        if str(sys.argv[1]).lower()[0] == "f":
            toFetch = False

        main(toFetch)

        print("""The GMTKN55 database has some inconsistent/missing data."""\
              """ Please see a generated file named 'Inconsistency'"""\
              """ and make manual changes to the Data.json file. """\
              """ The missing data can be fetched from the GMTKN55 database """\
              """ available at: https://www.chemie.uni-bonn.de/grimme/de/software/gmtkn/gmtkn55 """)
