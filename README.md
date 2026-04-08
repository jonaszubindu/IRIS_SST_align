# IRIS / SST Manual Alignment Workflow

This repository is now intentionally minimal. The workflow is:

1. Create a simple initial alignment.
2. Refine the SST pointing by hand in the interactive app.
3. Reuse that solved SST wideband WCS for narrowband products and overlap checks.

The codebase no longer keeps the earlier experimental affine, non-rigid, or multi-stage search scripts.

## Files

- [align_iris_sst.py](/Users/jonaszbinden/Documents/Playground/align_iris_sst.py)
  Resets the outputs to a clean initial alignment:
  - one global SST shift
  - optional IRIS shift
  - optional IRIS rotation

- [interactive_sst_manual_align.py](/Users/jonaszbinden/Documents/Playground/interactive_sst_manual_align.py)
  Plotly app for manual SST center placement on top of IRIS.

- [sunpy_alignment_check.py](/Users/jonaszbinden/Documents/Playground/sunpy_alignment_check.py)
  Rebuilds the simple SunPy blink GIF for the currently saved solution.

- [update_sst_nb_wcs.py](/Users/jonaszbinden/Documents/Playground/update_sst_nb_wcs.py)
  CLI wrapper that propagates the aligned SST WB spatial WCS into an NB cube.

- [plot_sst_iris_time_overlap.py](/Users/jonaszbinden/Documents/Playground/plot_sst_iris_time_overlap.py)
  Makes multi-timestep overlap plots for IRIS versus SST WB or SST NB wing images.

- [alignment_common.py](/Users/jonaszbinden/Documents/Playground/alignment_common.py)
  Shared FITS, timing, WCS, and plotting helpers.

## Inputs

- WB source: `/Users/jonaszbinden/Desktop/Align_IRIS_SST_proj/wb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_im.fits`
- NB source: `/Users/jonaszbinden/Desktop/Align_IRIS_SST_proj/nb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_cmapcorr_im.fits`
- IRIS source: `/Users/jonaszbinden/Desktop/Align_IRIS_SST_proj/iris_l2_20250619_072925_3660106834_SJI_2832_t000.fits`
- HMI reference: `/Users/jonaszbinden/Documents/Playground/data/hmi/hmi.sharp_720s.13354.20250619_083600_TAI.continuum.fits`

## Initial Alignment

Run:

```bash
python /Users/jonaszbinden/Documents/Playground/align_iris_sst.py --reset
```

By default this applies:

- SST shift: `(+3.60263, -3.59411)` arcsec
- IRIS shift: `(0.0, 0.0)` arcsec
- IRIS rotation: `0.0` deg

The outputs are:

- `/Users/jonaszbinden/Documents/Playground/outputs/wb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_im_aligned.fits`
- `/Users/jonaszbinden/Documents/Playground/outputs/iris_l2_20250619_072925_3660106834_SJI_2832_t000_aligned.fits`
- `/Users/jonaszbinden/Documents/Playground/outputs/alignment_report.json`

## Manual Refinement

Run:

```bash
python /Users/jonaszbinden/Documents/Playground/interactive_sst_manual_align.py
```

The app:

- opens in the browser
- shows IRIS and the projected SST WB image
- uses one mouse click to set the new SST center
- saves with `s` or the `Save` button
- enables `Quit And Close` after a successful save

Saving updates:

- `/Users/jonaszbinden/Documents/Playground/outputs/wb_3950_2025-06-19T08:35:11_08:35:11=0-76_corrected_im_aligned.fits`
- `/Users/jonaszbinden/Documents/Playground/outputs/sst_manual_adjustment.json`
- `/Users/jonaszbinden/Documents/Playground/outputs/alignment_report.json`
- `/Users/jonaszbinden/Documents/Playground/outputs/alignment_check_sunpy.gif`

## Narrowband WCS Propagation

The NB updater uses the aligned WB scan-by-scan corners and times, then builds the NB `WCS-TAB` coordinates as:

- `x, y` from the aligned WB cube
- `wavelength` from `WAVEMIN`, `WAVEMAX`, and `NAXIS3`
- `time` from the aligned WB scan times

Run:

```bash
python /Users/jonaszbinden/Documents/Playground/update_sst_nb_wcs.py
```

Important note:

The local NB file is still truncated, so writing a full aligned NB FITS here may fail until a complete cube is available. The propagation function is in place for a complete NB file.

## Time-Overlap Plots

To inspect the stability of the pointing over time:

```bash
python /Users/jonaszbinden/Documents/Playground/plot_sst_iris_time_overlap.py --sst-kind wb
```

For NB wing images on a complete cube:

```bash
python /Users/jonaszbinden/Documents/Playground/plot_sst_iris_time_overlap.py --sst-kind nb
```

Outputs:

- `/Users/jonaszbinden/Documents/Playground/outputs/sst_iris_time_overlap.png`
- `/Users/jonaszbinden/Documents/Playground/outputs/sst_iris_time_overlap_matches.json`

The NB mode uses the outermost wavelength planes as a proxy for the WB morphology.
