""" energy conversion factors """
__energy_conv_factors__ = {
    "Ry-Ha": 0.5,
    "Ha-Ry": 2.0,
    "Ha-eV": 27.211407952665,
    "eV-Ha": 0.03674929286053584,
    "Ry-eV": 13.6057039763325,
    "eV-Ry": 0.07349858572107168,
    "Ha-kcal": 627.5094741,
    "kcal-Ha": 0.0015936014375468055,
    "Ry-kcal": 313.75473705,
    "kcal-Ry": 0.003187202875093611,
    "eV-kcal": 23.06052943646173,
    "kcal-eV": 0.04336413883103952
}
__energyUnitStrsMap__ = {"r": "Ry", "h": "Ha", "a": "Ha", "e": "eV", "k": "kcal"}

""" length conversion factors """
__length_conv_factors__ = {
    "Angs-Bohr": 1.8897259885789,
    "Bohr-Angs": 0.5291772490000065,
}
__lengthUnitStrsMap__ = {"b": "Bohr", "au": "Bohr", "a": "Angs"}


def getEnergyUnitStr(unit):
    firstChar = unit.lower()[0]
    if firstChar not in list(__energyUnitStrsMap__.keys()):
        raise ValueError(
            """Invalid string for energy unit provided. """
            """Valid strings are: 'Rydberg' or 'Ry', 'Hartree' or 'Ha', """
            """ 'kcal/mol' or 'kcal', 'atomicUnit' or 'au', 'electronVolt' or 'eV'."""
        )
    return __energyUnitStrsMap__[firstChar]


def getLengthUnitStr(unit):
    shortStr = ""
    if len(unit) >= 2 and unit.lower()[0:2] == "au":
        shortStr = "au"
    else:
        shortStr = unit.lower()[0]

    if shortStr not in list(__lengthUnitStrsMap__.keys()):
        raise ValueError(
            """Invalid string for length unit provided. """
            """Valid strings are: 'Bohr' or 'b', 'Angstrom' or 'Angs', """
            """ 'atomicUnit' or 'au'."""
        )

    return __lengthUnitStrsMap__[shortStr]


def energyConversionFactor(inUnit, outUnit):
    inStr = getEnergyUnitStr(inUnit)
    outStr = getEnergyUnitStr(outUnit)
    if inStr == outStr:
        return 1.0

    else:
        return __energy_conv_factors__[inStr + "-" + outStr]


def lengthConversionFactor(inUnit, outUnit):
    inStr = getLengthUnitStr(inUnit)
    outStr = getLengthUnitStr(outUnit)
    if inStr == outStr:
        return 1.0
    else:
        return __length_conv_factors__[inStr + "-" + outStr]


def areSameEnergyUnits(a, b):
    return getEnergyUnitStr(a) == getEnergyUnitStr(b)


def areSameLengthUnits(a, b):
    return getLengthUnitStr(a) == getLengthUnitStr(b)
