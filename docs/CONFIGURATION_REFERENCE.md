# Configuration reference

The JSON configuration is the scientific model supplied to the common fitting
engine. Sample-specific Python builders record how revised JSON files were
constructed; they do not use a different fitting algorithm.

## Required top-level fields

- `frame_temperatures_C`: mapping from frame identifier to temperature.
- `profiles_to_compare`: any of `gaussian`, `lorentzian`, `pseudo_voigt`.
- `fit_settings`: background, multistart, resolution and detection settings.
- `windows`: list of independent local q-windows.

## Window fields

- `key` and `title`: stable identifier and display label.
- `q_min`, `q_max`: data interval used by the fit.
- `background_order`: normally 0 or 1.
- `profiles_to_compare`: optional per-window override.
- `eta_fixed`: optional locked pseudo-Voigt mixing value.
- `temperature_min_C`, `temperature_max_C`: optional applicability range.
- `peaks`: components included in the local model.

## Peak fields

- `key` and `label`: stable component identity.
- `q_guess`, `q_min`, `q_max`: initial centre and bounds in inverse angstrom.
- `fwhm_guess`, `fwhm_min`, `fwhm_max`: initial apparent FWHM and bounds.
- `anchor_frame` or `anchor_temperature_C`: optional reference for weak-feature
  stabilisation.
- `position_only`: centre/presence may be considered but width/area is
  intentionally suppressed.
- `reference_like_component`: nuisance/underlay component excluded from
  scientific result tables.
- `suppress_fwhm_interpretation`: retains the fit while explicitly preventing
  width interpretation.

Bounds are not measurement results. They encode the identity and physically
plausible range of the component and must be reviewed against the raw selected
frames for every new sample.
