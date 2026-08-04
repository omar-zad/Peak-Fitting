# Script map

## Fitting

- `fit_giwaxs_series.py`: the common fitting engine used for every sample and
  cut. It implements constrained model comparison, component detection,
  optional residual bootstrap, figures and numerical outputs.
- `run_ip_ccl_sensitivity.py`: a targeted Drop40 IP q approximately 0.306
  sensitivity workflow. It is not the general fitting entry point.

## Preparation and configuration

- `prepare_secondary_scan.py`: validates temperature, monitor and FR/IP/OOP
  alignment and writes prepared raw/normalised fitting tables.
- `build_drop40_nine_frame_inputs.py`: rebuilds the corrected nine-frame
  Drop40 inputs from external source data.
- `build_peakfit_cut_configs.py`: derives cut-specific Drop40 configs from a
  base definition.
- `build_in_situ_scan_configs.py`: constructs configs for an in-situ series.
- `build_secondary_sample_configs.py`: generic starting-config builder for a
  secondary scan; generated candidates require visual review.
- `build_scan587214_representative_configs.py`: representative Drop50 configs.
- `build_587217_candidate_v2.py`, `build_scan587225_candidate_v2.py`,
  `build_scan587225_candidate_v3.py`, `build_scan587241_candidate_v1.py`,
  `build_scan587250_candidate_v2.py` and `build_scan587250_candidate_v3.py`:
  audit trails for sample-specific revisions after visual diagnostics. These
  builders change the JSON model, not the fitting algorithm.

## Diagnostics

- `plot_key_frames.py`: unsmoothed selected-frame plots used before fitting to
  check candidate peaks and temperature regimes.

## Reporting

- `build_secondary_sample_report.py`: combines FR/IP/OOP tables and applies
  conservative numerical acceptance gates.
- `apply_manual_visual_qc.py`: applies a documented post-fit exclusion list.
- `build_nine_frame_qc_report.py`: combined Drop40 QC package.
- `build_ip_q0306_final_table.py`: targeted final Drop40 IP q approximately
  0.306 table.
- `summarize_drop40_cross_sector_sensitivity.py`: Drop40 FR/IP/OOP sensitivity
  summary.
- `build_secondary_summary_pdf.py`: optional PDF summary of secondary-sample
  reports.

## Calibration

- `sector_matched_calibrant_widths.py`: optional AgBh/LaB6 sector-matched
  instrumental-width diagnostic. It requires the packages in
  `requirements-calibration.txt` and external raw calibration files.
