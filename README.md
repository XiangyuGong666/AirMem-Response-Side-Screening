# AirMem Response-Side Screening

A standalone desktop application and Python workflow for the V19.0 response-side screening manuscript. It recomputes recurrent peak identification, CSD–MAC matching, equivalent added mass, and EFDD/RDT damping from the monitored wind and acceleration. The excitation-side software is a separate project.

## Data

From the [open monitoring dataset](https://doi.org/10.5281/zenodo.23994733), extract `AirMem_monitoring_10Hz_2022Q4.zip`. Choose the extracted `2022Q4` directory or its parent. A sibling `2023Q1` directory can stay in place; this workflow uses the manuscript period, October–December 2022. The program reads the 12-channel `analyacc_22_101112_dropbaseline_downsample10hz_12vertical.mat` and synchronized `analywindspeed_downsample_10hz.mat` without changing them. No manuscript image or precomputed measured-analysis output is distributed with the program.

The bare finite element model exports are bundled as scientific inputs for FE matching, geometry, signed modes and added-mass calculations. They are documented in [the model attribution note](src/airmem_response_screening/model_inputs/README.md). The program computes the modal mass quantities and corrected frequency library from these inputs; it does not launch an external FE solver.

## Windows download

For a ready-to-run copy, download [AirMem-Response-Side-Screening.exe] from [Release v0.1.2]. The single executable opens directly on 64-bit Windows and includes the Python runtime, desktop toolkit, and finite element inputs. Download the monitoring dataset separately.

[AirMem-Response-Side-Screening.exe]: https://github.com/XiangyuGong666/AirMem-Response-Side-Screening/releases/download/v0.1.2/AirMem-Response-Side-Screening.exe
[Release v0.1.2]: https://github.com/XiangyuGong666/AirMem-Response-Side-Screening/releases/tag/v0.1.2

## Run the analysis

Double-click the executable and follow these steps.

1. Choose **Select monitoring data** and select the local data directory.
2. Choose **Run monitoring calculation** and select a destination for the new results. The app creates a `Response-V19-...` subfolder and reports progress.
3. Select a figure and click **Generate selected figure**. The app displays it and saves an English PNG and PDF. To resume later, choose **Open completed run** and select the `Response-V19-...` results folder. Keep the matching acceleration file accessible for Figures 8, 13 and 20.

Use **Optional time-vector ZIP** only when supplying the original timestamps for Figure 3. The default is the public calendar index.

## Run from source

Use Python 3.11 or 3.12.

```powershell
python -m venv .venv
.venv/Scripts/python -m pip install .
.venv/Scripts/python -m airmem_response_screening
```

The CLI performs the same stages, with `plot` accepting one or more figure numbers (or `all` for scripted checks):

```powershell
.venv/Scripts/python -m airmem_response_screening calculate D:/Zenodo/2022Q4 results/v19-q4
.venv/Scripts/python -m airmem_response_screening plot results/v19-q4 7 12 21 22 23
```

## Figure coverage and verification

- Figures 3–23 can be generated. Figure 1 (methods flowchart) and Figure 2 (site photograph) are excluded. Figures 4–6 and 16 are calculated or drawn from the bundled FE inputs. The other figures are drawn from the monitoring run and its calculated intermediates.
- The public `paper1_calendar_2022Q4.csv` produces an explicitly labeled public-index version of Figure 3. The manuscript counts coverage from the original timestamps of 1,875 files. To reproduce that time basis, supply the original MAT directory via `--original-time`, or supply a ZIP of timestamp-named NPZ files containing the original `time` vectors via the GUI or CLI. The public-index version differed only on 2022-12-31 in the validation run, by 0.106444 h. The original time vectors are a separate input and are not bundled here.
- `run.json` records input and output SHA-256 checksums, the code fingerprints and method settings. `summary.json` contains group, mass and damping statistics. CSV files provide the full peak events, ranking, selected group, frequency comparisons (Table 3), quality decisions and wind bins (Table 4). `figure_manifest.json` records the figure input sources, their hashes, the completed run digest, and the Figure 3 time basis. An edited or missing intermediate is rejected rather than silently plotted. Run directories from before v0.1.2 must be recomputed.
- For the manuscript dataset, the expected audit values are 26,399 screened windows, 131,995 peaks and 118 groups. The selected group is ranked first at 0.550 Hz (1,261 peaks, 1,245 windows); calibrated added mass is 1.33334 kg/m². EFDD and RDT pass 1,233 and 653 windows, with medians of 0.550% and 0.752%; 648 windows pass both. These are comparison values only: they are not read by the analysis.

To compare a full run with the manuscript numerically, run `python scripts/check-v19.py results/v19-q4 --require-figures`. This check reads the outputs only.

The code is distributed under MIT. The [bundled FE input attribution](src/airmem_response_screening/model_inputs/README.md) and the [monitoring data record](https://doi.org/10.5281/zenodo.23994733) describe the separate terms for data and model inputs.

## Build the Windows executable

Use a clean Windows Python environment. Install the project and PyInstaller, then build:

```powershell
.venv/Scripts/python -m pip install . pyinstaller==6.20.0
powershell -ExecutionPolicy Bypass -File scripts/build-windows.ps1 -Python .venv/Scripts/python.exe -OutputRoot D:/AirMemBuild
```

The release artifact is `D:/AirMemBuild/dist/AirMem-Response-Side-Screening.exe`.
