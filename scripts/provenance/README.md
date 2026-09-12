# Provenance scripts

The scripts in this folder are kept as an audit trail, not as tools.

Each one records how a committed configuration or a dissertation table was
derived: which fitting windows were added, removed or re-bounded after a
visual review of the raw frames, and why. That reasoning is written in their
docstrings and comments, and the configurations they produced are committed
under `configs/`.

They read intermediate files from the private working tree (earlier candidate
configurations, prepared scan data, result folders) that are not part of this
repository, so they do not run from a fresh clone and are not covered by the
tests. They are deliberately left unmodified.

| Script | What it produced |
|---|---|
| `build_scan587214_representative_configs.py` | `configs/secondary_samples/scan_587214/representative/`, the seven-frame model-development configs that later served as the candidate-window template for the other secondary scans |
| `build_scan587214_comprehensive_v2.py` | the first comprehensive revision of scan 587214, superseded by v3 and kept for the audit trail |
| `build_scan587214_comprehensive_v3.py` | `configs/secondary_samples/scan_587214/comprehensive_v3/`, the current scan 587214 configurations |
| `build_587217_candidate_v2.py` | `configs/secondary_samples/scan_587217/` |
| `build_scan587225_candidate_v2.py`, `build_scan587225_candidate_v3.py` | `configs/secondary_samples/scan_587225/` (v3 imports v2 and applies the first-frame checkpoint revisions) |
| `build_scan587241_candidate_v1.py` | `configs/secondary_samples/scan_587241/` |
| `build_scan587250_candidate_v2.py`, `build_scan587250_candidate_v3.py` | `configs/secondary_samples/scan_587250/` |
| `build_ip_q0306_final_table.py` | the provisional Drop40 IP q ≈ 0.306 Å⁻¹ bootstrap table; it reads a specific 100-repeat bootstrap run |
| `build_secondary_summary_pdf.py` | the secondary-sample PDF summary; its per-sample notes record the visual-QC decisions taken for each scan |

The general-purpose tools that produced the remaining configurations
(`build_peakfit_cut_configs.py`, `build_secondary_sample_configs.py`,
`build_in_situ_scan_configs.py`) live in `scripts/preparation/` and are
tested.
