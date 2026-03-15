# File Descriptions

## sample_NNGGA_NNLDA.py
Header python code to be copied into the PySCF script. This header is supposed to contain the NNLDA and NNGGA related classes and functions.

## modules_env
File containing the module loads and the environment variables to be used inside the SLURM jobscript. 
Example
```
module load openmpi
module load cuda
module load cudnn
source ~/.bashrc
conda activate <my-conda-environment>
export <XYZ>=<ABC>
```

## GMTKN55.json
### Overview

**GMTKN55** (General Main-group Thermochemistry, Kinetics and Noncovalent Interactions, 55 subsets) is a comprehensive benchmark dataset for evaluating the accuracy of density functional approximations (DFAs) in quantum chemistry. It was developed by the Grimme group and covers a broad range of chemical properties across main-group elements.

This JSON file encodes the full GMTKN55 dataset including reference values and precomputed DFA results for **83 density functionals**, each evaluated with and without dispersion corrections.

**Reference publication:**
Goerigk, L.; Hansen, A.; Bauer, C.; Ehrlich, S.; Najibi, A.; Grimme, S. *Phys. Chem. Chem. Phys.* **2017**, *19*, 32184–32215.

---

### Top-Level Structure

The JSON file is a single dictionary with **55 keys**, one per benchmark subset. Each key is the subset name (e.g., `"W4-11"`, `"BH76"`, `"S66"`).

```json
{
  "W4-11": { ... },
  "G21EA": { ... },
  "BH76":  { ... },
  ...
}
```

**Total reactions across all subsets:** 1504

---

### Subset Dictionary

Each subset value is itself a dictionary keyed by **string-formatted integers** (`"1"`, `"2"`, ..., `"N"`), where N is the number of reactions in that subset.

```json
"W4-11": {
  "1": { ... },
  "2": { ... },
  ...
  "140": { ... }
}
```

---

### Reaction Entry

Each reaction entry is a dictionary with the following fields:

#### `Systems` — `list[str]`
A list of molecule/species identifiers (lowercase strings) involved in the reaction. These are short labels that correspond to the chemical systems used to define the reaction.

**Example:** `["h2", "h"]`

#### `Stoichiometry` — `list[int | float]`
A list of stoichiometric coefficients corresponding to each entry in `Systems`, in the same order. By convention:
- **Negative** values indicate **reactants**
- **Positive** values indicate **products**

The reaction energy is computed as the linear combination of system energies weighted by their stoichiometric coefficients.

**Example:** `[-1, 2]` means the reaction is: H₂ → 2 H

### `Ref.` — `float`
The high-accuracy reference (benchmark) value for the reaction energy or property, in **kcal/mol**. These values are used as the ground truth against which DFA results are compared.

**Example:** `109.493`

#### Functional keys — `dict`
All remaining keys are the names of density functionals (83 in total). Each maps to a sub-dictionary with **three dispersion correction variants**:

| Key | Description |
|---|---|
| `"without"` | Result with no dispersion correction |
| `"D3(0)"` | Result with Grimme's D3 dispersion correction using zero-damping |
| `"D3(BJ)"` | Result with D3 dispersion correction using Becke-Johnson damping |

Each variant value is either a `float` (in kcal/mol) or `null` if that combination was not available/computed.

**Example:**
```json
"PBE": {
  "without": 104.723,
  "D3(0)":   104.723,
  "D3(BJ)":  104.783
}
```

---

### Full Entry Example

```json
"W4-11": {
  "1": {
    "Systems":      ["h2", "h"],
    "Stoichiometry": [-1, 2],
    "Ref.":          109.493,
    "PBE": {
      "without": 104.723,
      "D3(0)":   104.723,
      "D3(BJ)":  104.783
    },
    "B3LYP": {
      "without": 110.423,
      "D3(0)":   110.423,
      "D3(BJ)":  110.513
    },
    ...
  }
}
```

This entry represents the atomization energy of H₂ into 2 H atoms. The reference value is 109.493 kcal/mol. PBE without dispersion predicts 104.723 kcal/mol (error: −4.770 kcal/mol).

---

### Units

All reaction energies and computed values are in **kcal/mol**.

---

### The 55 Subsets

The subsets span five main categories of chemical properties:

#### Thermochemistry (basic)
| Subset | Description |
|---|---|
| W4-11 | Atomization energies of small molecules |
| G21EA | Adiabatic electron affinities |
| G21IP | Adiabatic ionization potentials |
| DIPCS10 | Double ionization potentials of closed-shell systems |
| PA26 | Adiabatic proton affinities |
| SIE4x4 | Self-interaction error-related problems |
| ALKBDE10 | BDEs of alkaline and alkaline-earth metal compounds |
| YBDE18 | BDEs of transition metal complexes (Y-based) |
| AL2X6 | Dimerization energies of Al₂X₆ compounds |
| HEAVYSB11 | BDEs of heavy main-group compounds (Sb-based) |
| NBPRC | Oligomerization and H₂ activation for NH₃/BH₃ and PH₃/BH₃ |
| ALK8 | Dissociation of alkali metal compounds |
| RC21 | Reaction energies of 21 reactions with heavy elements |
| G2RC | Reaction energies of G2/97 test set |
| BH76RC | Reaction energies from BH76 barrier height set |
| FH51 | Reaction energies of 51 organic reactions |
| TAUT15 | Tautomerization energies |
| DC13 | Difficult cases (multireference, etc.) |
| MB16-43 | Dissociation of artificial molecules |

#### Thermochemistry (isomerization & reaction energies)
| Subset | Description |
|---|---|
| DARC | Diels-Alder reaction energies |
| RSE43 | Radical stabilization energies |
| BSR36 | Bond separation reactions of hydrocarbons |
| CDIE20 | Conformational/isomerization energies of dienes |
| ISO34 | Isomerization energies of small organic molecules |
| ISOL24 | Isomerization energies of large organic molecules |
| C60ISO | Isomerization energies of C₆₀ fullerene |
| PArel | Relative proton affinities |

#### Barrier Heights (kinetics)
| Subset | Description |
|---|---|
| BH76 | Hydrogen transfer, heavy atom transfer, nucleophilic substitution, etc. |
| BHPERI | Barrier heights of pericyclic reactions |
| BHDIV10 | Diverse reaction barrier heights |
| INV24 | Inversion barriers of amine and phosphine compounds |
| BHROT27 | Torsional/rotation barrier heights |
| PX13 | Proton exchange barriers |
| WCPT18 | Water-catalyzed proton transfer barriers |

#### Noncovalent Interactions
| Subset | Description |
|---|---|
| RG18 | Interaction energies of rare-gas dimers |
| ADIM6 | Interaction energies of alkane dimers |
| S22 | Noncovalent interactions in biological-sized molecules |
| S66 | Noncovalent interactions (diverse chemical space) |
| HEAVY28 | Noncovalent interactions of heavy element systems |
| WATER27 | Binding energies of water clusters |
| CARBHB12 | C–H···B hydrogen bonds |
| PNICO23 | Pnicogen bond interactions |
| HAL59 | Halogen bonding interactions |
| AHB21 | Anion–neutral hydrogen bonding |
| CHB6 | Cation–neutral hydrogen bonding |
| IL16 | Interaction energies in ionic liquid clusters |

#### Intramolecular Noncovalent Interactions & Conformational Energies
| Subset | Description |
|---|---|
| IDISP | Intramolecular dispersion interactions |
| ICONF | Intramolecular conformational energies |
| ACONF | Alkane conformational energies |
| Amino20x4 | Amino acid conformational energies |
| PCONF21 | Peptide conformational energies |
| MCONF | Melatonin conformational energies |
| SCONF | Sugar conformational energies |
| UPU23 | Uracil-containing RNA backbone conformational energies |
| BUT14DIOL | 1,4-butanediol conformational energies |

---

### The 83 Density Functionals

The functionals span all major rungs of Jacob's Ladder:

| Rung | Type | Examples |
|---|---|---|
| 1 | LDA | (not present) |
| 2 | GGA | PBE, BLYP, BP86, PW91, RPBE, ... |
| 3 | meta-GGA | TPSS, SCAN, M06L, MN15L, ... |
| 4 | Hybrid GGA | B3LYP, PBE0, M06, ... |
| 4 | Hybrid meta-GGA | M062X, TPSSh, BMK, MN15, ... |
| 4 | Range-separated hybrid | LC-ωhPBE, ωB97X-D3, ωB97X-V, HSE06, ... |
| 5 | Double hybrid | B2PLYP, B2GPPLYP, PWPB95, DSD-PBEP86, ... |

---

### Notes

- `null` values indicate that a particular functional/dispersion combination was not available in the original dataset.
- The `"D3(0)"` and `"D3(BJ)"` values for functionals that were parameterized with dispersion included by design (e.g., `B97-D3`, `wB97X-D3`) represent their canonical parametrization, while `"without"` may be `null` for those functionals.
- Stoichiometric convention (negative = reactants, positive = products) is consistent across all subsets.
- System name strings are lowercase shorthand labels; they do not follow a universal naming convention but are consistent within the dataset.

## BasicBarrierSystems.json 
   JSON file containing a sub-selection of GMTKN55 subsets. It contains the subset names and for each subset it contains the list of systems to consider. The format is
   ```json
   { "<subset1>": ["<subset1-sys1>", "<subset1-sys2>", ...], 
     "<susbet2>": ["<subset2-sys1>", "<subset2-sys2>", ...], 
     ...
    }
    ```
## Geom.json
    # GMTKN55 Geometry Dataset — README

### Overview

This JSON file contains the **molecular geometries** for all chemical systems appearing across the 55 subsets of the GMTKN55 benchmark dataset. Each molecule is stored with its charge, spin multiplicity, and Cartesian atomic coordinates.

- **Coordinate units:** Ångström (Å)
- **Total subsets:** 55
- **Total molecule entries:** 2548 (some molecules appear in more than one subset)
- **Largest molecule:** 81 atoms (`i4e` in `ISOL24`)

---

### Top-Level Structure

The JSON file is a single dictionary with **55 keys**, one per GMTKN55 subset. The subset names match exactly those used in the companion `GMTKN55.json` energetics file (e.g., `"W4-11"`, `"BH76"`, `"S66"`).

```json
{
  "W4-11":  { ... },
  "G21EA":  { ... },
  "BH76":   { ... },
  ...
}
```

---

### Subset Dictionary

Each subset value is a dictionary keyed by **molecule name** (lowercase string label). These labels correspond directly to the entries in the `"Systems"` lists of each reaction in `GMTKN55.json`.

```json
"W4-11": {
  "h2":           { ... },
  "h":            { ... },
  "acetaldehyde": { ... },
  ...
}
```

The number of molecules per subset varies (e.g., 152 in `W4-11`, 198 in `S66`). Note that **355 molecule names appear in more than one subset** — in those cases, the same name may hold the same or a closely related geometry depending on context (e.g., `"al"` appears in both `W4-11` and `G21IP`).

---

### Molecule Entry

Each molecule entry is a dictionary with three fields:

#### `charge` — `int`
The total electric charge of the molecule. Possible values in this dataset: `-1`, `0`, `1`, `2`.

**Examples:**
- Neutral molecule: `0`
- Anion: `-1` (e.g., `EA_10` in `G21EA`)
- Cation: `1` (e.g., `Ap` in `PA26`)
- Dication: `2` (e.g., `be_2+` in `DIPCS10`)

#### `mult` — `int`
The spin multiplicity of the molecule, defined as $2S + 1$ where $S$ is the total spin. Possible values in this dataset: `1`, `2`, `3`, `4`.

| Value | Spin state | Example |
|---|---|---|
| `1` | Singlet | Most closed-shell molecules |
| `2` | Doublet | Radicals, open-shell atoms (e.g., H atom) |
| `3` | Triplet | Biradicals (e.g., O atom, NH) |
| `4` | Quartet | High-spin open-shell species |

#### `geom` — `list[list[str]]`
A list of atoms, where each atom is represented as a list of four strings:

```
[element_symbol, x, y, z]
```

| Element | Type | Description |
|---|---|---|
| `element_symbol` | `str` | Standard atomic symbol (e.g., `"H"`, `"C"`, `"Al"`) |
| `x` | `str` | x-coordinate in Ångström, formatted as a decimal string |
| `y` | `str` | y-coordinate in Ångström, formatted as a decimal string |
| `z` | `str` | z-coordinate in Ångström, formatted as a decimal string |

> **Note:** Coordinates are stored as **strings**, not floats. Convert to `float` before numerical use (e.g., `float(x)`).

---

### Full Entry Examples

**Diatomic molecule (H₂, neutral singlet):**
```json
"h2": {
  "charge": 0,
  "mult": 1,
  "geom": [
    ["H", "0.000000", "0.000000",  "0.370946"],
    ["H", "0.000000", "0.000000", "-0.370946"]
  ]
}
```

**Atom (H radical, neutral doublet):**
```json
"h": {
  "charge": 0,
  "mult": 2,
  "geom": [
    ["H", "0.000000", "0.000000", "0.000000"]
  ]
}
```

**Polyatomic molecule (acetaldehyde, neutral singlet, 7 atoms):**
```json
"acetaldehyde": {
  "charge": 0,
  "mult": 1,
  "geom": [
    ["C",  "0.000000",  "0.463403",  "0.000000"],
    ["O",  "1.205540",  "0.374658",  "0.000000"],
    ["H", "-0.484415",  "1.458495",  "0.000000"],
    ["C", "-0.936133", "-0.711830",  "0.000000"],
    ["H", "-0.376435", "-1.644095",  "0.000000"],
    ["H", "-1.583337", "-0.660551",  "0.878440"],
    ["H", "-1.583337", "-0.660551", "-0.878440"]
  ]
}
```

**Dication (Be²⁺, singlet):**
```json
"be_2+": {
  "charge": 2,
  "mult": 1,
  "geom": [
    ["Be", "0.000000", "0.000000", "0.000000"]
  ]
}
```

---

### Relationship to GMTKN55.json

This file is the geometry companion to `GMTKN55.json`. The two files are linked as follows:

- Each reaction in `GMTKN55.json` has a `"Systems"` list of molecule name strings (e.g., `["h2", "h"]`) and a `"Stoichiometry"` list of signed coefficients.
- The same molecule name strings are the keys within the matching subset in this geometry file.
- To fully characterise a reaction, look up each system name in the corresponding subset of `Geom.json` to retrieve its structure, charge, and multiplicity.

**Example lookup for W4-11 reaction 1:**

```
GMTKN55.json → W4-11 → "1" → Systems: ["h2", "h"], Stoichiometry: [-1, 2]
Geom.json    → W4-11 → "h2"  : charge=0, mult=1, 2 atoms
Geom.json    → W4-11 → "h"   : charge=0, mult=2, 1 atom
```

---

### Notes

- Coordinates are **zero-centered** per molecule but no specific orientation convention is enforced globally.
- Single atoms are included as one-element `geom` lists (coordinates `[0.0, 0.0, 0.0]`).
- All coordinate strings use six decimal places of precision.
- Molecule name strings are lowercase and may include characters such as `+`, `-`, or digits (e.g., `"h2+"`, `"be_2+"`, `"EA_c-"`).
