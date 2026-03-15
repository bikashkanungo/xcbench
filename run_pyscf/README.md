# Description of different scripts and how to run them

## CreatePyscfScript.py

### Overview

`CreatePyscfScript.py` automates the generation and submission of PySCF DFT calculations for a user-defined set of molecular systems. Given a JSON input file, it:

1. Selects a range of molecules from a systems JSON and a geometry JSON.
2. Generates a self-contained PySCF Python script for each molecule, assembling it from a user-supplied header, a molecule-build block, an optional SCF guess block, and the main SCF block.
3. Writes a single SLURM job script that allocates all required nodes at once and launches one `srun` step per molecule in parallel (a multi-step SLURM job), then immediately submits it with `sbatch`.

---

### Usage

```bash
python CreatePyscfScript.py inp.json
python CreatePyscfScript.py inp.json START END
```

| Argument | Required | Description |
|---|---|---|
| `inp.json` | Yes | Path to the JSON input file (must contain `.json`) |
| `START` | No | Integer index overriding `sysRange.start` in the JSON |
| `END` | No | Integer index overriding `sysRange.end` in the JSON |

Pass `-h` or `--help` to print usage information.

---

### What the Script Produces

#### Per-molecule PySCF scripts

For each molecule, a directory `<rootDir>/<subset>/<sysName>/` is created and a Python script named `<sysName><sysFilePostfix>.py` is written there. The script is assembled in up to three parts:

**1. Header** (`header` file, always included)
A user-supplied Python file prepended verbatim. Typically contains imports (`pyscf`, `numpy`, custom XC classes like `NNGGA`/`NNLDA`, etc.) and any environment setup.

**2. Molecule build block** (always included)
Defines the PySCF `Mole` object and instantiates the KS object:
- Sets `mol.atom` from the geometry data (Cartesian coordinates in Ångström)
- Sets `mol.charge`, `mol.spin` (= multiplicity − 1), and `mol.basis`
- Uses `dft.RKS` for closed-shell (spin=0) and `dft.UKS` for open-shell systems

**3. Pre-SCF guess block** (optional, when `useGuess.flag` is `true`)
Runs a cheaper SCF first to generate a density matrix as the initial guess for the main SCF. The guess functional is always mapped to either PBE (for GGA-type XC strings) or SPW92 (for LDA-type). The SCF solver can be `diis` or `newton`. After convergence, the density matrix `dm` is passed to the main SCF via `mfl.init_guess = dm`.

**4. Main SCF block** (always included)
Configures and runs the target SCF calculation:
- Sets the XC functional (see [Supported XC Functionals](#supported-xc-functionals))
- Sets `max_cycle = 100` and the convergence tolerance
- Applies the occupancy scheme if specified (see [Occupancy](#occupancy))
- Runs `mfl.kernel()` and saves the density matrix to a text file:
  - Closed-shell: `DM_<sysFilePostfix>`
  - Open-shell: `DM0_<sysFilePostfix>` (alpha) and `DM1_<sysFilePostfix>` (beta)

#### SLURM job script

A single SLURM batch script is written to `<jobscript>_<start>_<end>` (suffixed with `_2`, `_3`, ... if the file already exists). It:
- Requests `nodesPerSys × N` nodes in total, where N is the number of systems
- Sources the `modulesAndEnv` file for module loads and environment variables
- Sets `OMP_NUM_THREADS`
- For each molecule, `cd`s to its directory and launches an independent `srun` step with `&`, so all systems run simultaneously within the same allocation
- Ends with a `wait` to hold the job until all steps finish
- Copies the `ptcFile` (neural network parameter file) to each system's directory before launching

The script is submitted immediately via `os.system("sbatch " + jobscript)`.

---

### Input JSON Reference

Below is a description of every field in the input JSON, illustrated with the example `inpBasicBarrier.json`.

#### File and Directory Paths

| Key | Type | Description |
|---|---|---|
| `rootDir` | `str` | Root directory under which all per-system subdirectories are created. Resolved to an absolute path. |
| `sysJSON` | `str` | Path to the systems JSON file. This is a subset-structured dictionary (same format as `GMTKN55.json`) listing the molecular systems to be computed. |
| `geomJSON` | `str` | Path to the geometry JSON file (`Geom.json`), containing the charge, multiplicity, and Cartesian coordinates for each system. |
| `header` | `str` | Path to a plain Python file prepended verbatim to every generated PySCF script. |
| `modulesAndEnv` | `str` | Path to a plain text/shell file containing `module load` commands and environment variable exports, sourced into the SLURM job script. |
| `jobscript` | `str` | Base filename for the generated SLURM job script. The actual filename will be `<jobscript>_<start>_<end>`. |
| `outFile` | `str` | Filename used to redirect the stdout/stderr of each `srun` step (i.e., `&> <outFile>`). |
| `stdOut` | `str` or `null` | If a string, redirects the Python-level stdout (system selection report) to this file. If `null`, output goes to stdout. |

#### PySCF Calculation Parameters

| Key | Type | Description |
|---|---|---|
| `basis` | `str` | Basis set string passed to `mol.basis` (e.g., `"def2-qzvp"`). |
| `convTol` | `float` | SCF convergence tolerance for the main SCF (`mfl.conv_tol`). |
| `printDM` | `bool` | Whether to save the density matrix. Currently wired as `true`; the density matrix is always written. |
| `sysFilePostfix` | `str` | String appended to the molecule name to form the PySCF script filename and the density matrix output filename. |

#### XC Functional (`xc`)

| Sub-key | Type | Description |
|---|---|---|
| `type` | `str` | XC functional identifier. See [Supported XC Functionals](#supported-xc-functionals). |
| `ptcPath` | `str` | Path to the neural network parameter file (`.ptc`), required for `NNLDA` and `NNGGA`. Copied to each system directory at runtime. Ignored for standard functionals. |
| `tol` | `float` | Numerical tolerance passed to the `NNLDA`/`NNGGA` constructor. |
| `sthres` | `float` | Density gradient threshold (`sthres`) passed to the `NNGGA` constructor. |

##### Supported XC Functionals

| `type` value | Description |
|---|---|
| `"PBE"` | Standard PBE GGA functional via PySCF's built-in `'PBE,PBE'` |
| `"PW92"` | PW92 LDA functional via PySCF's `'SPW92'` |
| `"SCAN"` | SCAN meta-GGA; uses `'PBE,PBE'` as the pre-SCF guess XC |
| `"R2SCAN"` | r²SCAN meta-GGA via PySCF's `'R2SCAN'` |
| `"NNLDA"` | Neural-network LDA functional; requires `ptcPath`. Loaded via a custom `NNLDA` class and `eval_xc_lda` hook. |
| `"NNGGA"` | Neural-network GGA functional; requires `ptcPath`. Loaded via a custom `NNGGA` class and `eval_xc_gga` hook. |

#### SCF Guess (`useGuess`)

Controls whether a cheaper pre-SCF calculation is run first to generate a density matrix initial guess for the main SCF.

| Sub-key | Type | Description |
|---|---|---|
| `flag` | `bool` | If `true`, the pre-SCF guess step is included in the script. |
| `xcStr` | `str` | XC functional for the guess SCF. GGA-type strings (`"PBE"`, `"NNGGA"`, `"SCAN"`, `"R2SCAN"`) map to `'PBE,PBE'`; LDA-type (`"PW92"`, `"NNLDA"`) map to `'SPW92'`. |
| `method` | `str` | SCF solver for the guess step: `"diis"` (default DIIS) or `"newton"` (second-order Newton–Raphson). |
| `convTol` | `float` | Convergence tolerance for the guess SCF (typically looser than the main SCF). |

#### Occupancy (`occupancy`)

Controls the orbital occupancy scheme applied to the main SCF.

| Sub-key | Type | Description |
|---|---|---|
| `type` | `str` or `null` | `"frac"` for fractional occupancy (`scf.addons.frac_occ`), `"fermi"` for Fermi smearing (`scf.addons.smearing_`), or `null` to use the default integer occupancy. |
| `T` | `float` | Electronic temperature in Kelvin, used only when `type` is `"fermi"`. Converted to Hartree via $\sigma = k_B T$ (with $k_B = 3.1668 \times 10^{-6}$ Ha/K). |

#### System Range (`sysRange`)

Selects a contiguous slice of molecules from the flattened list of all systems in `sysJSON`. Systems are enumerated in the order they appear in the JSON (subset by subset, then molecule by molecule within each subset).

| Sub-key | Type | Description |
|---|---|---|
| `start` | `int` | Index of the first system to include (inclusive). |
| `end` | `int` | Index of the last system to include (exclusive). |

These values can be overridden from the command line: `python CreatePyscfScript.py inp.json START END`.

#### SLURM Parameters (`slurm`)

| Sub-key | Type | Description |
|---|---|---|
| `jobname` | `str` | Base SLURM job name; the actual name is `<jobname>_<start>_<end>`. |
| `account` | `str` | SLURM account (`#SBATCH -A`). |
| `partition` | `str` or `""` | SLURM partition (`#SBATCH -p`). Omitted if empty or `null`. |
| `queue` | `str` or `""` | SLURM QOS/queue (`#SBATCH -q`). Omitted if empty or `null`. |
| `t` | `str` | Wall-clock time limit (e.g., `"4:00:00"`). |
| `nodesPerSys` | `int` | Number of nodes allocated per system. Total nodes requested = `nodesPerSys × N`. |
| `tasksPerNode` | `int` | Passed to `srun` as `--ntasks-per-node`. |
| `cpusPerTask` | `int` | Passed to `srun` as `--cpus-per-task`. |
| `threads` | `int` | Value exported as `OMP_NUM_THREADS`. |
| `extras` | `list[str]` | Additional `#SBATCH` directives, written verbatim (e.g., `["-C cpu"]`). |
| `email` | `str` | Email address for `BEGIN`/`END` notifications. |

---

### Output Directory Structure

```
<rootDir>/
└── <subset>/
    └── <sysName>/
        ├── <sysName><sysFilePostfix>.py   # PySCF script
        ├── <ptcFile>                       # Copied NN parameter file
        ├── DM_<sysFilePostfix>             # Density matrix (closed-shell)
        ├── DM0_<sysFilePostfix>            # Alpha density matrix (open-shell)
        └── DM1_<sysFilePostfix>            # Beta density matrix (open-shell)
```

The SLURM job script is written to the directory from which the script is invoked:
```
<jobscript>_<start>_<end>
```

---

### Example

Using `inpBasicBarrier.json`:

```bash
python CreatePyscfScript.py inpBasicBarrier.json
```

This will:
- Load systems 0–1099 from `BasicBarrierSystems.json` and their geometries from `Geom.json`
- For each system, generate a PySCF script using `def2-qzvp`, `NNGGA` as the target functional with parameters from `NNGGA_UEG.ptc`, a PBE/Newton pre-SCF guess, fractional occupancy, and a convergence tolerance of 1e-7
- Create a SLURM job requesting `1 × N` nodes on the `regular` queue on the `m3059` account, with one `srun` per system running in parallel, and submit it

To override the system range from the command line:

```bash
python CreatePyscfScript.py inpBasicBarrier.json 0 200
```

---

### Internal Functions

| Function | Description |
|---|---|
| `getSysData(sysJSON, geomJSON, start, end)` | Loads the `[start:end]` slice of systems from the systems JSON, looks up each system's geometry data in the geometry JSON, and returns a `dict` keyed by `"subset:sysName"`. |
| `getMolBuildFile(sysData, basisStr, outFname)` | Writes the `mol = gto.Mole()` block and KS instantiation to a temporary file. |
| `getPreSCFFile(sysData, xcStr, method, convTol, outFname)` | Writes the pre-SCF guess block to a temporary file. |
| `getSCFFile(sysData, xcStr, convTol, occupancy, ptcFile, tol, sthres, printDM, postfix, outFname)` | Writes the main SCF block (XC setup, occupancy, kernel call, DM save) to a temporary file. |
| `createScripts(systems, ...)` | Orchestrates directory creation and script assembly for all systems by concatenating the header and temporary files. |
| `run(inp, stdoutF)` | Top-level function: parses the input dict, calls `getSysData`, `createScripts`, writes the SLURM job script, and submits it. |




## EnergyDiff.py — README

### Overview

`EnergyDiff.py` reads PySCF DFT output files for a set of molecular systems, computes reaction energies, and evaluates error metrics against GMTKN55 reference data. It reports per-subset mean errors (ME) and mean absolute errors (MAE) for a user-specified list of DFAs (drawn from the precomputed values in `GMTKN55.json`) and for the NN functional whose energies are parsed from PySCF output files. It also computes two weighted thermochemical mean absolute deviation scores (WTMAD-1 and WTMAD-2). Results are printed to stdout and two files are written: `weights.json` and `errMat`.

---

### Usage

```bash
python EnergyDiff.py inp.json
```

The script must be run from the directory that contains the per-system output subdirectories, i.e., the same `rootDir` used when generating calculations with `CreatePyscfScript.py`. The expected directory structure is:

```
<cwd>/
└── <subset>/
    └── <sysName>/
        └── <outFileToParse>
```

---

### Input JSON Reference

| Key | Type | Description |
|---|---|---|
| `sysJSON` | `str` | Path to a systems JSON file (e.g., `BasicBarrierSystems.json`). Each subset entry must contain an `"IDs"` list specifying which reaction indices to include. This allows working with a "diet" subset of the full GMTKN55. |
| `allDataJSON` | `str` | Path to `GMTKN55.json`, containing reference energies and precomputed DFA results for all reactions. |
| `outFileToParse` | `str` | Filename of the PySCF output file within each system's directory. Must match `outFile` used in `CreatePyscfScript.py`. |
| `subsets` | `list[str]` | List of GMTKN55 subset names to include in the analysis. |
| `DFAs` | `list[str]` | List of DFA names (must match keys in `GMTKN55.json`) to compare against the NN functional. |

---

### What the Script Does

#### 1. System and reaction selection

`getSysNames` and `getReactionsData` build the working set of molecules and reactions by cross-referencing `sysJSON` (which specifies the reaction IDs of interest per subset) with `allDataJSON` (which holds the full GMTKN55 data). Only the reactions listed under each subset's `"IDs"` key in `sysJSON` are used.

#### 2. Energy parsing (`parseEnergy`)

For each molecule, the PySCF output file at `<cwd>/<subset>/<sysName>/<outFileToParse>` is read. Two convergence strategies are tried in order:

**Primary — converged SCF line:** The file is scanned in reverse for a line containing `"converged SCF energy"`. The energy (in Hartree) is parsed from the text after `=` and converted to kcal/mol using `HARTREE_TO_KCAL_PER_MOL = 627.5096080305927`.

**Fallback — near-convergence RMSE check:** If no converged line is found and `rmseTol` is not `None`, the last `NLast` (default: 10) SCF cycle lines containing `"cycle="` are parsed for the `E=` field. If the RMSE of these energies falls below `rmseTol` (default: 5×10⁻⁴ Ha), the molecule is treated as effectively converged and the mean energy is used.

Systems that pass neither check are absent from the returned `energies` dict. Any subset containing at least one such unconverged system is automatically added to `subsetsToSkip` and excluded from the error analysis.

#### 3. Error computation

For each reaction in each considered subset, the NN reaction energy is computed as:

```
ΔE_NN = Σ_i stoichiometry[i] × energy[system_i]
```

The signed error and absolute error relative to the GMTKN55 reference value `Ref.` are accumulated into per-subset ME and MAE for both the NN and each DFA in `DFAs`. DFA reaction energies are read from the `"without"` dispersion variant in `GMTKN55.json`.

#### 4. Weighting schemes

Two WTMAD variants are computed, following the conventions of the original GMTKN55 paper. Both depend on `avgRef[subset]`, the mean of `|Ref.|` over all reactions in a subset.

**WTMAD-1:** Each subset's weight is inversely proportional to its average reference magnitude, with clamping to prevent extreme subsets from dominating:

| Condition on `avgRef` (kcal/mol) | Weight |
|---|---|
| > 75 | 0.1 / N_subsets |
| < 7.5 | 10.0 / N_subsets |
| otherwise | 1.0 / N_subsets |

**WTMAD-2:** Each subset's weight is proportional to the number of reactions it contains and inversely proportional to its average reference magnitude, normalized by the grand average reference energy across all reactions:

```
w2[subset] = (allAvgRef / nallReactions) × (nReactions[subset] / avgRef[subset])
```

where `allAvgRef` is the arithmetic mean of `avgRef[subset]` over all considered subsets.

#### 5. Output

**Stdout:**
- List of subsets being considered
- Per-subset line: MAE for every method (NN + all DFAs), WTMAD-2 weight, reaction count, and average reference magnitude
- WTMAD-1 and WTMAD-2 totals for all methods (WTMAD-0 is printed but always zero — see Issues)
- Per-DFA comparison: subsets where NN is worse than / better than each reference DFA (currently hardcoded to PBE and SCAN — see Issues)

**`weights.json`:** WTMAD-1 (`w1`) and WTMAD-2 (`w2`) weights for each considered subset.

**`errMat`:** A plain-text matrix of MAE values with shape `(len(DFAs)+1) × (N_subsets+2)`. Rows are DFAs + NN; columns are the considered subsets followed by WTMAD-1 and WTMAD-2.

---

### Internal Functions

| Function | Description |
|---|---|
| `parseEnergy(sysNames, rootoutfilename, rmseTol, NLast)` | Parses PySCF output files; returns a dict mapping `"subset:sysName"` to energy in kcal/mol for converged systems only. |
| `getSysNames(dietReactions, allData, subsets)` | Returns the unique list of `"subset:sysName"` strings needed across all selected reactions. |
| `getReactionsData(dietReactions, allData, subsets)` | Returns a nested dict of reaction data for the selected subset/reaction-ID combinations. |
| `main(dietReactions, allData, subsets, outFileToParse, DFAs)` | Orchestrates energy parsing, error computation, weighting, and all output. |

---
