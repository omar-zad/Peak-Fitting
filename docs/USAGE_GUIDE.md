# Usage guide

## 1. Understand the inputs

The fitting input is a table with one q column and one intensity column per
frame. A JSON configuration maps frame identifiers to temperature and defines
each local fitting window.

The fitter operates on already integrated 1D line cuts. The `.poni`, mask and
raw detector files are therefore not required to reproduce an ordinary fit.
They remain relevant to the upstream pyFAI/pygix integration and separate
instrumental-resolution work.

All commands below write under `local_results/`, which Git ignores.

## 2. Inspect profiles before fitting

For a secondary scan prepared with `prepare_secondary_scan.py`, first plot the
selected unsmoothed frames:

```bash
python scripts/diagnostics/plot_key_frames.py \
  --processed-dir local_results/prepared/<sample>/scan_000000 \
  --scan-id 000000 \
  --sample-label "Sample label" \
  --frames 0,10,20 \
  --outdir local_results/000000/key_frame_plots
```

The processed folder must contain `scan000000_<cut>_<N>frames_norm.txt` for
each cut and `scan000000_frame_manifest.csv`, which is what the preparation
script writes.

Use this plot to choose local windows, component counts and initial bounds.
Manual peak labels are initial hypotheses, not final reported positions.

## 3. First-frame checkpoint

Use `--stage first` for one diagnostic frame. It intentionally stops after one
frame and writes `FIRST_FRAME_REVIEW.txt` and `first_frame_diagnostic.png`.

Check that:

1. the background follows the local baseline;
2. the total fit follows the visible feature;
3. components correspond to plausible peaks or shoulders;
4. residuals do not contain repeated peak-shaped structure;
5. centres and widths are not repeatedly at bounds;
6. neighbouring components retain their identities.

## 4. Selected-frame screen

Use `--stage selected --frames ... --bootstrap 0` after the checkpoint. The
profile comparison is repeated across the active frames in each window unless
the configuration explicitly locks the profile family.

Run one cut at a time: FR, inspect, then IP, inspect, then OOP.

## 5. Profile selection

Gaussian, Lorentzian and area-normalised pseudo-Voigt profiles are fitted to
the same q-window with the same component count, bounds and background order.
BIC is combined across all active frames in that window. Lowest BIC wins,
except that a simpler family is retained when pseudo-Voigt's improvement is
less than 2 BIC units.

Pseudo-Voigt is a Gaussian/Lorentzian mixture. Its eta value is shared across
the active frames for that window. Profile selection is therefore not a
per-frame majority vote.

## 6. Peak detection and reporting

Each component is tested by removing it and refitting. The default evidence
levels combine leave-one-component-out delta BIC with integrated-area SNR:

| Level | Delta BIC | Area SNR |
|---|---:|---:|
| tentative | >= 2 | >= 2 |
| detected | >= 6 | >= 3 |
| strong | >= 10 | >= 5 |

The detection status is assigned from the preliminary all-free fit and
re-derived from the final fit only in the anchor frame, so a component labelled
tentative in another frame can carry a final area SNR above 3. The CSV shows
the preliminary delta BIC and the final SNR side by side.

Shape-bound components are unresolved or rejected. Weak components may be
stabilised using a reference-frame shape, but those borrowed q/FWHM values are
not independent measurements. A component whose FWHM hits its lower bound in
the preliminary fit falls back to a borrowed shape unless it is declared
`position_only`, in which case it is fitted freely and only its position may be
reported.

## 7. Visual QC and reporting

Generate the combined FR/IP/OOP report:

```bash
python scripts/reporting/build_secondary_sample_report.py \
  --scan 000000 \
  --sample-label "Sample label" \
  --manifest local_results/prepared/<sample>/scan_000000/scan000000_frame_manifest.csv \
  --results-root local_results/000000/key_frame_fits \
  --outdir local_results/000000/key_frame_report
```

`--results-root` must contain `FR/`, `IP/` and `OOP/` run folders. The report
writes `accepted_peak_positions.csv`, `accepted_peak_areas.csv`,
`tentative_or_excluded_peaks.csv`, `all_peak_fits_with_gates.csv`,
`selected_profile_families.csv`, `REPORT.md` and `report_provenance.json`.

If a numerically accepted position is visibly unreliable, document it in a
manual exclusion CSV with the columns `cut, frame, window, peak, reason` and
run:

```bash
python scripts/reporting/apply_manual_visual_qc.py \
  --report-dir local_results/000000/key_frame_report \
  --exclusions path/to/manual_qc_exclusions.csv
```

The preferred final table is
`accepted_peak_positions_after_visual_qc.csv` when generated.

## 8. All-frame and bootstrap runs

Only use `--stage all --confirm-first-fit` when a continuous transition or
kinetic conclusion requires every frame. Use `--bootstrap 0` for configuration
development. Reserve bootstrap repeats for a small number of final conclusions
whose uncertainty materially matters.

For long series the fitter writes a representative overview figure plus a
`<figure>_pages/` folder holding every frame, so nothing is hidden by the
figure size.
