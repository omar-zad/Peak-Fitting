# Usage guide

## 1. Understand the inputs

The fitting input is a table with one q column and one intensity column per
frame. A JSON configuration maps frame identifiers to temperature and defines
each local fitting window.

The fitter operates on already integrated 1D line cuts. The `.poni`, mask and
raw detector files are therefore not required to reproduce an ordinary fit.
They remain relevant to the upstream pyFAI/pygix integration and separate
instrumental-resolution work.

## 2. Inspect profiles before fitting

For a secondary scan, first plot the selected unsmoothed frames:

```bash
python scripts/diagnostics/plot_key_frames.py \
  --processed-dir path/to/processed \
  --scan-id 000000 \
  --sample-label "Sample label" \
  --frames 0,10,20 \
  --outdir local_results/000000/key_frame_plots
```

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

Shape-bound components are unresolved or rejected. Weak components may be
stabilised using a reference-frame shape, but those borrowed q/FWHM values are
not independent measurements.

## 7. Visual QC and reporting

Generate the combined FR/IP/OOP report:

```bash
python scripts/reporting/build_secondary_sample_report.py \
  --scan 000000 \
  --sample-label "Sample label" \
  --manifest path/to/frame_manifest.csv \
  --results-root local_results/000000/key_frame_fits \
  --outdir local_results/000000/key_frame_report
```

If a numerically accepted position is visibly unreliable, document it in a
manual exclusion CSV and run:

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
