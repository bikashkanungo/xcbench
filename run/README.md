#Explanation of key-value pairs in inp.json

1. `"mode"`: Defined the mode of operation. Valid types are`"create"`, `"run"`, `"both"`.\ 

    1. "create": Creates each of the system folders and the neccessary files within in 
       (parametes.prm, coordinates.inp, domainVectors.inp, MLXC.inp, MLXC.ptc, etc)
    2. "run": Launches DFT-FE for each of the systems
    3. "both": Does both "create" and "run" in that order.

2. `"rootDir"`: The parent directory within which all the system directories are created or assumed to exist.

3. `"inpParamsFile"`: Speficies an input DFT-FE parameters which is used as a base parameters file 
    to create the system specific parameters file 

4. `"outParamsFile"` Name to be given to the final system specific DFT-FE parameters file.\ 
   Example: "parameters.prm"
  
5. `"outCoordsFile"`: Name to be given to the system specific coordinates file.\
   Example: "coordinates.inp"
  
6. `"outDomainFile"`: Name to be give to the system specific domain vectors file.\
   Example: "domainVectors.inp"

7. `"outFile"`: Name of the file to redirect stdout and stderr while running DFT-FE.\
   Example: "out"

7. `"vaccuumBuffer"`: Width of the vacuum to be considered to specify the domain vectors.
   Recommended to use 25.0.
  
8. `"xc"` : Dictionary of the following form specifying the exchange-correlation functional and its related parameters. 
   ```python
   {"type": "<xc-type>", "modelPath": "<path for MLXC model.inp file>", "ptcPath": "<path for MLXC .ptc file>"}
   ```
   The `xc-type` (e.g., GGA-PBE, MLXC-NNGGA, etc.) are listed in DFT-FE manual.\
   If not using any MLXC, then the `"modelPath"` and `"ptcPath"` values can be `null`.

9. `"dispersion"`: Specifies the dispersion type. Valid types are: `"D3"`, `"D4"`, and `null`.

10. `"D3Damping"`: Specifies the damping type for D3 dispersion type. Valid values: "0" for no damping,  "BJ" for Becke-Johnson damping. It is used only when `"dispersion": "D3"`.
  
11. `"sysRanges"`: Dictionary of the form
    ```python
    {"start": <start-id>, "end": <end-id>}
    ``` 
    where `start-id` and `end-id` specify the range of systems to consider. The indices start from 0. The start-id is included but the `end-id` is not.If `start-id` being `null` is assumed to be 0 and `end-id` being `null` is assumed to be till the end of all systems provided via `"sysJSON"`. This option allows one to run DFT-FE on a selection of systems durind a batch run.

12. `"dftfeExec"`: Specifies the path to the DFT-FE exectuable. 

13. `"srunParams"`: String to specify the parameters to the `srun` launcher. It can be `null` to let the SLURM decide based on the batch script.\
    Example: `"-np 128 --exclusive"` to run it on 128 processors and with `"exclusive"` option. 

12. "sysJSON": JSON file containing the list of systems to consider. The assumed format is a dictionary is as follows:
    ```python
      { 
        "subset-name-1": ["system-name-1", system-name-2",..."system-name-2"],
        "subset-name-2": ["system-name-1", system-name-2",..."system-name-2"],
        ...
        so on and so forth for other subsets
        ...
        }
    ```
13. "geomJSON": JSON file containing the geometries of all the systems across all the subsets. The assumed format is a nested dictionary of the form:
    ```python
      {
        "subset-name-1": {
          "system-name-1": {
            "charge": <charge-value>,
              "mult": <spin-multiplicity>,
              "geom": [
                [
                  "<atomic-symbol-1>",
                  "<x-coordinate>",
                  "<y-coordinate>",
                  "<z-coordinate>"
                ],
                [
                  "<atomic-symbol-2>",
                  "<x-coordinate>",
                  "<y-coordinate>",
                  "<z-coordinate>"
                ],
                ...
                ...
                [
                  "<atomic-symbol-N>",
                  "<x-coordinate>",
                  "<y-coordinate>",
                  "<z-coordinate>"
                ]
              ]
          },
          "system-name-2": {
            "charge": <charge-value>,
            "mult": <spin-multiplicity>,
            "geom": [
              [
                "<atomic-symbol-1>",
                "<x-coordinate>",
                "<y-coordinate>",
                "<z-coordinate>"
              ],
              [
                "<atomic-symbol-2>",
                "<x-coordinate>",
                "<y-coordinate>",
                "<z-coordinate>"
              ],
              ...
              ...
              [
                "<atomic-symbol-N>",
                "<x-coordinate>",
                "<y-coordinate>",
                "<z-coordinate>"
              ]
            ]
          },
          ...
          so and so forth for other systems in subset-name-1
          ...
        },
        {
          ....
          subset-name-2 data
          ....
        },
        ...
        so on and so forth for other subsets
        ...
      }
    ```
