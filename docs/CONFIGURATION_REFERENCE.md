# Configuration reference

The JSON configuration is the scientific model supplied to the common fitting
engine. The builders in `scripts/preparation/` and the audit-trail scripts in
`scripts/provenance/` record how each committed file was constructed; they do
not use a different fitting algorithm.

## Required top-level fields

- `frame_temperatures_C`: mapping from frame identifier to temperature.
- `profiles_to_compare`: any of `gaussian`, `lorentzian`, `pseudo_voigt`.
- `fit_settings`: background, multistart, resolution and detection settings.
- `windows`: list of independent local q-windows.

Optional top-level fields such as `analysis_variant`, `analysis_scope`,
`representative_frames` and `series_metadata` are informational: they record
what a configuration is for, which frames were reviewed and where it came from,
and do not change the fit. Pass the frames to fit explicitly with `--frames`.

## Window fields

- `key` and `title`: stable identifier and display label.
- `q_min`, `q_max`: data interval used by the fit.
- `background_order`: normally 0 or 1.
- `profiles_to_compare`: optional per-window override.
- `eta_fixed`: optional locked pseudo-Voigt mixing value.
- `temperature_min_C`, `temperature_max_C`: optional applicability range. A
  frame outside it is skipped for that window.
- `frame_min`, `frame_max`: optional zero-based, inclusive acquisition-frame
  range that is applied in addition to the temperature range. Use it only for
  isothermal holds whose pattern changes with time rather than temperature,
  for example the two halves of the scan 587214 46 °C hold. A window without
  these keys behaves exactly as before they were introduced.
- `candidate_window_requires_visual_review`: informational flag written by the
  secondary-scan builder; it marks windows that were transferred from another
  sample and not yet reviewed against this sample's own frames.
- `peaks`: components included in the local model.

## Peak fields

- `key` and `label`: stable component identity.
- `q_guess`, `q_min`, `q_max`: initial centre and bounds in inverse angstrom.
- `fwhm_guess`, `fwhm_min`, `fwhm_max`: initial apparent FWHM and bounds.
- `anchor_frame` or `anchor_temperature_C`: optional reference for weak-feature
  stabilisation.
- `position_only`: the component is fitted freely and may support presence and
  q position, but its width and area are intentionally not reported. Use it for
  lines only one to three q-steps wide, which would otherwise fall back to a
  borrowed shape when their FWHM hits the lower bound.
- `reference_like_component`: nuisance/underlay component excluded from
  scientific result tables.
- `suppress_fwhm_interpretation`: retains the fit while explicitly preventing
  width interpretation.
- `provisional_assignment`: informational flag marking a component whose
  assignment has not been confirmed.

Bounds are not measurement results. They encode the identity and physically
plausible range of the component and must be reviewed against the raw selected
frames for every new sample.

## Committed configurations

### `configs/drop40/`

All fifteen files are written by `scripts/preparation/build_peakfit_cut_configs.py`
from `base/peakfit_config.json`, and the test suite checks that they still are.

| File | Role |
|---|---|
| `base/peakfit_config.json` | shared starting model for all three cuts |
| `base/peakfit_config_OOP.json` | OOP one-broad-component low-q intermediate |
| `peakfit_config_FR_revised.json`, `peakfit_config_IP_revised.json`, `peakfit_config_OOP_revised.json` | the retained nine-frame models |
| `sensitivity/peakfit_config_FR_2component.json`, `sensitivity/peakfit_config_FR_3component_candidate.json` | FR low-q component-count alternatives |
| `sensitivity/peakfit_config_IP_1component.json`, `sensitivity/peakfit_config_IP_3component.json` | IP low-q component-count alternatives (README comparison) |
| `sensitivity/peakfit_config_OOP_lowq_quadratic_sensitivity.json` | OOP low-q quadratic-background check |
| `sensitivity/peakfit_config_{FR,IP,OOP}_pipi_highT_{one,two}.json` | high-temperature pi-pi component-count checks |

Two further files in `sensitivity/` were written by hand rather than by the
builder: `peakfit_config_IP_q0480_cross_sector.json` and
`peakfit_config_OOP_q0459_cross_sector.json`, the reciprocal cross-sector test
summarised by `summarize_drop40_cross_sector_sensitivity.py`.

### `configs/secondary_samples/`

One folder per scan holding the reviewed configuration for each cut. The
`series_metadata` block of each file records the key frames, the anchor frame
and, where a builder was used, the template it started from. Paths that point
into the private working tree are written as `<private-working-tree>/...`.

| Scan | Files | Origin |
|---|---|---|
| 587178, 587184, 587207 | `scan<id>_<cut>_candidate_keyframes.json` | built by `build_secondary_sample_configs.py`, then revised by hand after visual review; the `candidate_revision` block records the changes |
| 587214 | `representative/` and `comprehensive_v3/` | `representative/` is the seven-frame window library used as the template for the other scans; `comprehensive_v3/` is the current scan-specific analysis, which uses `frame_min`/`frame_max` gating and `position_only` lines |
| 587217, 587225, 587241, 587250 | `scan<id>_<cut>_candidate_keyframes.json` | written by the scan-specific builders in `scripts/provenance/` |

`scan_registry.example.json` is the template for the registry that
`prepare_secondary_scan.py` reads. The real registry contains local paths and
is Git-ignored.
