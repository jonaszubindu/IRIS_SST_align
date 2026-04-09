# IRIS / SST Alignment Workflow

This repository aligns SST wideband and narrowband products with IRIS SJI data, using HMI as the reference frame for the global coordinate system.

The workflow is intentionally simple:

1. Build a coarse initial alignment.
2. Refine IRIS manually against HMI.
3. Refine SST manually against IRIS.
4. Regenerate WB comparison products.
5. Optionally build the aligned NB cube if enough disk space is available.

## Repository Layout

- `align_iris_sst.py`
  Creates the initial aligned IRIS and SST WB products.

- `interactive_iris_hmi_align.py`
  Manual Plotly/Dash app for refining IRIS against HMI.

- `interactive_sst_manual_align.py`
  Manual Plotly/Dash app for refining SST against IRIS.

- `plot_vs_hmi.py`
  Creates HMI comparison figures for IRIS or SST.

- `plot_sst_iris_time_overlap.py`
  Creates time-overlap checks between SST and IRIS.

- `sunpy_alignment_check.py`
  Builds the SunPy blink GIF alignment check.

- `update_sst_nb_wcs.py`
  Propagates the aligned SST WB WCS into the NB cube.

- `alignment_common.py`
  Shared paths, FITS helpers, WCS logic, and plotting helpers.

- `run_full_pipeline.sh`
  Runs the full workflow, including the interactive steps.

## Installation

Install from the repository root:

```bash
pip install -e .
```

This installs command-line entry points such as:

- `iris-sst-align-init`
- `iris-sst-align-fetch-hmi`
- `iris-sst-align-run`
- `iris-sst-align-iris-app`
- `iris-sst-align-sst-app`
- `iris-sst-align-plot-hmi`
- `iris-sst-align-plot-overlap`
- `iris-sst-align-blink`
- `iris-sst-align-update-nb`

## Configuration On Another Machine

You can configure the repo without editing the code by setting environment variables:

- `IRIS_SST_ALIGN_DATA_DIR`
  Directory containing the SST, IRIS, and helper files.

- `IRIS_SST_ALIGN_HMI_PATH`
  Optional explicit path to the local HMI FITS file.

- `IRIS_SST_ALIGN_JSOC_EMAIL`
  Registered JSOC email address used when the HMI file needs to be downloaded automatically through `drms`.

If you prefer editing the code directly, the relevant paths are in `alignment_common.py`.

These source-file paths are derived from `DATA_DIR`, so in most cases changing `IRIS_SST_ALIGN_DATA_DIR` is enough:

- `SST_WB_SOURCE_PATH`
- `SST_NB_SOURCE_PATH`
- `IRIS_SOURCE_PATH`
- `STIC_UTILS_PATH`

You normally do not need to change:

- `WORKDIR`
  It is derived automatically from the location of `alignment_common.py`.

- `OUTPUT_DIR`
  It is derived automatically as `WORKDIR / "outputs"`.

## Required Input Files

The code expects these source files:

- `wb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_im.fits`
- `nb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_cmapcorr_im.fits`
- `iris_l2_20250619_072925_3660106834_SJI_2832_t000.fits`
- `STIC_pyfun_utils.py`

The HMI reference file is expected at:

- `data/hmi/hmi.sharp_720s.13354.20250619_083600_TAI.continuum.fits`

If it is missing, the code will try to download it automatically from JSOC via `drms` using `IRIS_SST_ALIGN_JSOC_EMAIL` or `JSOC_EMAIL`.

## Running The Workflow

Run from the repository root:

```bash
cd /path/to/IRIS_SST_align
bash run_full_pipeline.sh
```

After installation, you can also run:

```bash
iris-sst-align-run
```

The script will:

1. Check the main input files.
2. Remove old output products.
3. Run the initial coarse alignment.
4. Launch the interactive IRIS app.
5. Launch the interactive SST app.
6. Regenerate WB comparison products.
7. Ask whether you want to continue, rerun interactive steps, or exit.
8. Optionally generate NB comparison plots.
9. Optionally build the aligned NB cube.

## Manual Commands

If you want to run steps separately:

Initial alignment:

```bash
iris-sst-align-init --reset
```

IRIS manual alignment:

```bash
iris-sst-align-iris-app
```

SST manual alignment:

```bash
iris-sst-align-sst-app
```

WB comparison plots:

```bash
iris-sst-align-plot-hmi --sst-kind wb --num-samples 3
iris-sst-align-plot-overlap --sst-kind wb --num-samples 3
iris-sst-align-blink
```

NB comparison plots:

```bash
iris-sst-align-plot-hmi --sst-kind nb --num-samples 3
iris-sst-align-plot-overlap --sst-kind nb --num-samples 3
```

NB aligned cube:

```bash
iris-sst-align-update-nb
```

Explicit HMI download:

```bash
iris-sst-align-fetch-hmi --email you@example.org
```

## Outputs

The main outputs are written into `outputs/`.

Typical files are:

- `outputs/wb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_im_aligned.fits`
- `outputs/iris_l2_20250619_072925_3660106834_SJI_2832_t000_aligned.fits`
- `outputs/alignment_report.json`
- `outputs/compare_to_hmi_wb.png`
- `outputs/sst_iris_time_overlap_wb.png`
- `outputs/alignment_check_sunpy.gif`

If NB products are enabled and disk space allows:

- `outputs/compare_to_hmi_nb.png`
- `outputs/sst_iris_time_overlap_nb.png`
- `outputs/nb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_cmapcorr_im_aligned.fits`

## Disk Space Note

The aligned NB cube is large. Building it can require substantially more free disk space than the final file size alone, because FITS updates may temporarily need additional write room during flush/resize.

If disk space is tight, it is safer to:

- skip the NB aligned cube locally
- generate only WB and NB comparison plots
- or build the NB aligned cube on a larger remote machine
