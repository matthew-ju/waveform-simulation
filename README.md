# DLOPy: Horizontal Sensor Orientation (BK network)

Python script to estimate **horizontal seismometer orientation** by cross-correlating observed surface waves with theoretical predictions. Supports single-station runs or batch runs over many stations, and writes per-station figures plus a **summary CSV**. Easily configurable parameters through YAML file.

---

## Setup

```bash
# clone
git clone git@github.com:matthew-ju/waveform-simulation.git
cd waveform-simulation

# conda env (Python 3.11) + deps
conda create -n bsl python=3.11 -y
conda activate bsl
conda install -c conda-forge obspy numpy scipy matplotlib pandas pyyaml tqdm -y
```

---

## Configuration
1. Download and run the first helper script to produce a fresh `stations.txt` (one station per line):
```
python parse.py
```
Example `stations.txt`:
```
ADAM
ALVW
YUBA
```
2. Download the R.gv dispersion tables:
```
R.gv.10.txt  R.gv.15.txt  R.gv.20.txt  R.gv.25.txt
R.gv.30.txt  R.gv.35.txt  R.gv.40.txt
```

3. Download and configure `DLOPy.yml`
### Key configurable parameters at a glance

| Key | Meaning | Typical value |
| :-- | :------ | :------------ |
| `multiple` | Batch run over `sta_file` vs. single-station run via `sta` | `true` for batch |
| `net`, `sta`, `sta_file` | Network and station selection | `BK`, `"YUBA"`, `"stations.txt"` |
| `cha`, `com`, `loc` | Channel mask, component mask, location code | `BH?`, `BH?`, `00` |
| `Rdir`, `outdir` | Where dispersion tables live and where outputs go | `./DLOPy`, `./out` |
| `cat_client`, `wf_client` | Catalog and waveform data sources | `IRIS` / `NCEDC` |
| `time1`, `time2` | Analysis time window | `"YYYY-MM-DD HH:MM:SS"` |
| `minmag`, `mindeg_sw`, `maxdeg_sw`, `maxdep_sw` | Event QC | `6.5`, `5.0`, `175.0`, `150` |
| `constsave`, `finplot`, `savecat` | Save behavior and figures | `true` |
   

---


## How to Run
Run from within the window folder that contains `DLOPy.yml`, `stations.txt`, and `Rdir/`:

```bash
conda activate bsl
python3 DLOpy.py
```

Run secondary helper script `visualizations.py` utilizing the outputted `summary.csv`.
```bash
python3 visualizations.py
```

---
## Interpreting results
### 3 plots are generated for each station:
1. Orientation-Correlation (Scatter Plot)
   `[station_ID]_v2_cluster.pdf`
  * Raw DLOPy result, correlating its horizontal orientation in degrees (azimuth) with the cross-correlation value
  * Solid red line shows its mean azimuth from true north, with `0ª` being most optimal
  * Dotted red lines show the bootstrap confidence interval (CI) from the mean
2. Spread of Azimuth (Histogram)
    `[station_ID]_boostrap_means_hist.pdf`
  * Distribution of the initial spread of azimuth
  * Solid red line shows the mean spread of azimuth from true north, with `0º` being most optimal
  * Dotted red lines show the 95% CI for this mean
3. Bootstrapped Mean Distribution (Histogram)
   `[station_ID]_boostrap_hist.pdf`
  * Distribution of the mean azimuth calculated from 10,000 bootstrap resamples
  * The distribution should be approximately bell-shaped, showing the stability and precision of the estimated mean azimuth from the bootstrap procedure
  * Solid red line shows the mean spread of azimuth from true north, with `0º` being most optimal
  * Dotted red lines show the 95% CI for this mean

> Note: All calculations are used only with data points above the 95% confidence interval

### Dataset summarizing results `summary.csv` with typical columns:
| Field | Description| Example |
| :-- | :------ | :------------ |
| Network | Seismic network code | BK |
| Station | Station ID | ADAM |
| Location | Location code | 00 |
| Mean Azimuth | Mean orientation estimate in degrees relative to true north | 358.3 |
| CI95 Low Interval | Lower bound of 95% CI for azimuth | 355.7 |
| CI95 High Interval | Upper bound of 95% CI for azimuth | 0.9 |
| Data Points Used | Number of measurements/boostraps used in estimate | 256 |
| Within5? | Flag whether the mean azimuth is within 5º of the target/reference | Yes/No |
| AbsOffset | Absolute difference (in degrees) between mean azimuth and true north | 1.7 |

* Small `AbsOffset` within 5º in confidence interval -> likely aligned closer to 
* Large `AbsOffset` not within 5º not in confidence interval -> likely needs rotation ~`Mean Azimuth`

### 3 additional summarizing plots are generated for each overall run
1. Station Pass/Fail (Bar Chart)
   * High-level overview of the overall success rate of the stations in a network
   * Compares the number of stations for which the mean azimuth was within 5º of true north (pass) and those for the mean azimuth exceeded 5º (fail)
   * Grouped bar chart, where each bar shows the total number AND proportion of stations that "passed" and "failed"
2. Distribution of Errors in Degrees (Histogram & Box Plot)
  * Statistical distribution of absolute offset in degrees
  * All stations are binned to visualize distribution and key statistics like median, quartiles, and outliers
3. Failed Stations Details (Bar Chart)
   * Detailed analysis of the "failing" stations (only stations with mean azimuths exceeding 5º were considered)
   * Ranks the stations in order of deferring from true north, in degrees

---
## Re-running for a new time window
Create a new window subfolder, copy configurations and tables, adjust YAML file, and re-run:
```bash
# example: new window for 2024-01 to 2026-01
mkdir -p trial11/2024-01_to_2026_01
cp DLOPy.yml trial11/2024-01_to_2026_01/
cp -r DLOPy trial11/2024-01_to_2026_01/    # contains R.gv.*.txt
cp stations.txt trial11/2024-01_to_2026_01/

# edit trial11/2024-01_to_2026_01/DLOPy.yml:
#   - time1, time2
#   - outdir (point it to this window’s output folder)
cd trial11/2024-01_to_2026_01
python3 parse.py      # optional: refresh stations.txt
python3 DLOpy.py
```

---
## Troubleshooting
* `stations.txt` not found
  * Ensure `sta_file: stations.txt` exsists in the current working directory; re-generate with `python parse.py`
* `R.gv.<xx>.txt` not found
  * Confirm all `R.gv.10 … R.gv.40` files are present and `Rdir` points to that folder.
* Imports flagged yellow/`ModuleNotFoundError`
  * Use the Conda Python and select it in your editor:
  ```bash
  conda activate bsl
  which python
  python3 -c "import obspy, numpy, scipy, matplotlib, pandas, yaml"
  ```
* Service/network issues (script is unable to pull from IRIS online server)
  * Retry later.

---
## Citation & License
This workflow adapts ideas and supporting files from [BSL TOOLKIT (DLOPy)](https://github.com/sylvster/BSL_TOOLKIT) -- please credit the original project and respect its license.
Data services used: [IRIS](https://service.iris.edu ), [NCEDC](https://ncedc.org/web-services-home.html)
