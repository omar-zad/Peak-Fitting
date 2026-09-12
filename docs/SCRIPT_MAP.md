# Script map

Every script below runs from a fresh clone with the packages in
`requirements.txt`, unless stated otherwise. All of them write under the
Git-ignored `local_results/` folder by default.

## Fitting

- `fit_giwaxs_series.py`: the common fitting engine used for every sample and
  cut. It implements constrained model comparison, component detection,
  optional residual bootstrap, figures and numerical outputs. `--config` and
  `--outdir` are always required.
- `run_ip_ccl_sensitivity.py`: the Drop40 IP q ≈ 0.306 Å⁻¹ sensitivity audit.
  It drives the fitter through six background, window and profile-family
  variants of the retained IP configuration and writes a self-describing
  package to `local_results/IP_CCL_sensitivity/`. It is not a general entry
  point.

## Preparation and configuration

- `prepare_secondary_scan.py`: validates temperature, monitor and FR/IP/OOP
  alignment for one scan listed in `configs/secondary_samples/scan_registry.json`
  and writes the raw and d5i-normalised fitting tables, a frame manifest, a
  key-frame suggestion list and a provenance file. Acquisition-format
  specific; see the README's Scope section.
- `build_drop40_nine_frame_inputs.py`: rebuilds the corrected nine-frame
  Drop40 inputs in `example_data/drop40/` from the external source data
  (`docs/NINE_FRAME_DATA_CORRECTION.md`). Requires the raw exports, which are
  not in the repository.
- `build_peakfit_cut_configs.py`: derives every configuration under
  `configs/drop40/` from `configs/drop40/base/peakfit_config.json`: the
  retained FR, IP and OOP configurations and the sensitivity variants. The test
  suite checks that it regenerates the committed files exactly.
- `build_in_situ_scan_configs.py`: adapts the retained Drop40 windows to a new
  in-situ series by replacing only the frame-temperature map and the anchor
  temperatures.
- `build_secondary_sample_configs.py`: builds candidate configurations for a
  secondary scan from the scan 587214 representative window library, or, with
  `--freeze-from-results`, locks each window to the profile family selected in
  a completed key-frame run.

## Diagnostics

- `plot_key_frames.py`: unsmoothed selected-frame plots used before fitting to
  check candidate peaks and temperature regimes.

## Reporting

- `build_nine_frame_qc_report.py`: combines the nine-frame FR, IP and OOP runs,
  applies the reporting gates and writes the Drop40 QC package.
- `compare_component_models.py`: checks that saved candidate runs use the same
  q-range and point count, then exports a BIC comparison table and figure.
- `summarize_drop40_cross_sector_sensitivity.py`: summarises the reciprocal
  q ≈ 0.459/0.480 Å⁻¹ IP/OOP cross-sector test against the nine-frame baseline
  runs.
- `build_secondary_sample_report.py`: combines FR/IP/OOP key-frame tables for a
  secondary scan and applies the conservative numerical acceptance gates.
- `apply_manual_visual_qc.py`: applies a documented post-fit exclusion list
  (columns `cut, frame, window, peak, reason`) to the accepted positions.

## Calibration

- `sector_matched_calibrant_widths.py`: optional AgBh/LaB6 sector-matched
  instrumental-width diagnostic. It requires the packages in
  `requirements-calibration.txt` and external raw calibration files.

## Provenance

`scripts/provenance/` holds the one-off scripts that record how the committed
scan-specific configurations and two dissertation tables were derived. They
read intermediate files that are not in the repository, so they are kept as an
audit trail rather than as tools; its README lists what each one produced.
