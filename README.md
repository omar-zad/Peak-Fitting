# Constrained GIWAXS peak fitting

This private research repository contains the 1D GIWAXS line-profile fitting
workflow used for a dissertation project. It supports staged fitting of FR,
IP and OOP profiles using constrained Gaussian, Lorentzian and pseudo-Voigt
models, visual residual checks and conservative result reporting.

The repository is intentionally private. It contains small processed example
line cuts, but no raw detector images, Nexus/HDF5 files or unrestricted source
data.

## What the workflow does

For each configured q-window, the fitter:

1. fits all configured peak components plus a local polynomial background;
2. compares Gaussian, Lorentzian and pseudo-Voigt families using BIC across
   the active frames in that window;
3. favours a simpler profile when pseudo-Voigt improves BIC by less than 2;
4. tests each component by comparing the full model with a reduced model;
5. reports peak position, d-spacing, integrated area and apparent FWHM when
   the corresponding numerical and resolution checks pass;
6. saves overlays, individual components, residuals and an auditable table of
   all accepted, tentative and rejected fits.

No smoothing, vertical intensity offsets, MCMC or unconstrained exhaustive
search is used. The fit is a bounded nonlinear least-squares fit with a small
number of deterministic multistarts.

Profile selection is made per fitting window using all active frames. It is
not a majority vote over independently selected profiles for each peak.

## Installation

Python 3.10 or later is required. Python 3.11 or 3.12 is recommended for a
straightforward installation.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The calibrant-resolution script has additional packages:

```bash
python -m pip install -r requirements-calibration.txt
```

Those optional packages are not required for ordinary fitting of prepared 1D
q-intensity profiles.

## Reproduce the included Drop40 checkpoint

Run the deliberately short first-frame checkpoint:

```bash
python scripts/fitting/fit_giwaxs_series.py \
  example_data/drop40/drop40_FR_9frames_norm.txt \
  --config configs/drop40/peakfit_config_FR_revised.json \
  --outdir local_results/drop40_FR_checkpoint \
  --stage first \
  --bootstrap 0
```

Review these files before continuing:

- `local_results/drop40_FR_checkpoint/figures/first_frame_diagnostic.png`
- `local_results/drop40_FR_checkpoint/FIRST_FRAME_REVIEW.txt`
- `local_results/drop40_FR_checkpoint/figures/fits_<window>.png`
- `local_results/drop40_FR_checkpoint/fit_quality.csv`

A verified reference run is committed under
`example_results/drop40_FR_checkpoint/` for comparison. New runs go under the
Git-ignored `local_results/` directory so machine-specific metadata is not
accidentally committed.

Only after the diagnostic has been accepted, run all nine frames:

```bash
python scripts/fitting/fit_giwaxs_series.py \
  example_data/drop40/drop40_FR_9frames_norm.txt \
  --config configs/drop40/peakfit_config_FR_revised.json \
  --outdir local_results/drop40_FR_all \
  --stage all \
  --confirm-first-fit \
  --bootstrap 0
```

Replace `FR` with `IP` or `OOP` in both the input and configuration filename
to run the other cuts.

## Run selected frames from another prepared scan

```bash
python scripts/fitting/fit_giwaxs_series.py \
  path/to/scan587225_FR_135frames_norm.txt \
  --config configs/secondary_samples/scan_587225/scan587225_FR_candidate_keyframes.json \
  --outdir local_results/scan_587225/FR \
  --stage selected \
  --frames 0,10,19,21,24,27,30,40,70,134 \
  --bootstrap 0
```

Run FR first, inspect its overlays and residuals, then repeat for IP and OOP.
Do not assume a configuration from one sample is automatically valid for a
different sample.

## Prepare a new secondary scan

Copy the registry template and replace its placeholder paths with local paths:

```bash
cp configs/secondary_samples/scan_registry.example.json \
  configs/secondary_samples/scan_registry.json
```

Then run:

```bash
python scripts/preparation/prepare_secondary_scan.py --scan 000000
```

This validates temperature and frame alignment, matching q-grids, FR/IP/OOP
column completeness and positive d5i monitor values. It writes raw and
d5i-normalised fitting tables plus a provenance manifest. Raw experimental
data remain outside this repository.

## Important outputs

| Output | Meaning |
|---|---|
| `peak_parameters.csv` | Complete component table, including fitted and reportable values |
| `fit_quality.csv` | Optimiser, residual, BIC and bound diagnostics |
| `model_selection.csv` | Gaussian/Lorentzian/pseudo-Voigt comparison by window |
| `fitted_curves.csv` | Numerical data, total fit, background and components |
| `figures/fits_<window>.png` | Raw data, total fit, components and residuals |
| `config_used.json` | Exact configuration snapshot used by the run |
| `run_metadata.json` | Input, frames, seed and software versions |
| `run_summary.txt` | Concise profile-selection and quality summary |

For secondary-sample reports, use `accepted_peak_positions_after_visual_qc.csv`
when present. A nonblank internal `q0_fit_Ainv` is not automatically an
accepted scientific result; the conservative gate uses `q0_reported_Ainv`.

## Scientific limits

- Confidence intervals are conditional on the selected peak/background model.
- The current profiles do not contain propagated pointwise measurement errors;
  empirical local noise is used for comparison and diagnostics.
- Apparent FWHM is not instrument-deconvolved.
- Do not convert secondary-sample FWHM values to CCL, crystallite size or
  strain without an appropriate instrumental-resolution treatment.
- A radial 1D profile alone cannot establish face-on/edge-on orientation or
  mosaicity; retain the corresponding 2D/azimuthal evidence.

Further details are in [docs/USAGE_GUIDE.md](docs/USAGE_GUIDE.md),
[docs/METHODS_SUMMARY.md](docs/METHODS_SUMMARY.md),
[docs/CONFIGURATION_REFERENCE.md](docs/CONFIGURATION_REFERENCE.md) and
[docs/SCRIPT_MAP.md](docs/SCRIPT_MAP.md).

## Private sharing

This repository is for invite-only dissertation collaboration. See
[SHARING_NOTICE.md](SHARING_NOTICE.md). Do not change the repository visibility
or redistribute data without approval from the researcher, supervisor and any
applicable institution or beamline data policy.
