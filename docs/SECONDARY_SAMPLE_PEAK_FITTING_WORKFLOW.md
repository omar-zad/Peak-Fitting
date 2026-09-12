# Secondary-sample GIWAXS peak-fitting workflow

## Aim and scientific boundary

Use this workflow for the non-Drop40 samples to obtain reproducible peak
presence, fitted peak position (`q0`), d-spacing (`d = 2*pi/q0`), integrated
area, IP/OOP sector comparison, and changes with time/temperature.

Do not use this workflow to claim CCL, crystallite size, strain or disorder for
the secondary samples. FWHM remains a fitting diagnostic only.

No smoothing, vertical offset or MCMC is used. All initial runs use
`--bootstrap 0`.

## Status of each scan

The reviewed configuration for every scan below is committed under
`configs/secondary_samples/scan_<id>/`. The fitted results and the
dissertation workbooks derived from them are held in the author's private
working tree, not in this repository.

- Scan 586875 is excluded from the current dissertation peak-fitting workflow
  by decision of the author.
- Scan 587214 (drop50) has completed the eight selected-frame screen (0, 25,
  50, 57, 60, 71, 85, 99) in FR, IP and OOP with the scan-specific
  `comprehensive_v3` configurations and `--bootstrap 0`, after a full visual QC
  of the v2 overlays. The earlier v1 package is superseded: its six windows
  were tuned around frame 57 and hid many visible features. v3 reports 204
  accepted positions plus presence-only and tentative tiers. Frame 0 is
  amorphous.
- Scan 587178 has completed the eight confirmed key-frame screen in FR, IP and
  OOP using sample-specific candidate configs, `--bootstrap 0`, visual fit QC,
  and a conservative report. It is not an all-frame or CCL/FWHM analysis.
- Scan 587184 has completed the eight confirmed key-frame screen in FR, IP and
  OOP with sample-specific candidate-v2 configs, `--bootstrap 0`, full visual
  QC and an auditable manual exclusion.
- Scan 587207 has completed the nine confirmed key-frame screen in FR, IP and
  OOP with temperature-regime-specific candidate-v5 configs, `--bootstrap 0`,
  full visual QC and auditable manual exclusions for unstable FR high-q
  component centres. Frame 45 is labelled final/slight cooling, not a
  separately established cooled state.
- Scan 587217 has completed the ten confirmed key-frame screen in FR, IP and
  OOP with candidate-v2 configs and `--bootstrap 0`. One background-sensitive
  OOP shoulder was removed by visual QC, leaving 89 dissertation-facing q
  positions. Frame 35 is late cooling at 73.9 °C.
- Scan 587225 has completed the ten confirmed key-frame screen in FR, IP and
  OOP with revised candidate-v3 configs and `--bootstrap 0`. Reference-like
  underlays are excluded by construction and three dubious positions were
  removed by visual QC, leaving 135 dissertation-facing q positions.
- Scan 587241 has completed the thirteen confirmed key-frame screen in FR, IP
  and OOP with candidate-v1 configs and `--bootstrap 0`. The transition occurs
  between the selected 48.8 and 61.6 °C frames; 150 positions passed the
  conservative gate, with asymmetric or under-resolved lines retained only as
  cautious centre/presence observations.
- Scan 587250 has completed the eight confirmed key-frame screen in FR, IP and
  OOP with candidate-v3 configs and `--bootstrap 0`. Strong asymmetric high-q
  features are retained for provisional centre comparison only, not width
  interpretation.

## Confirmed key frames

| Scan | Condition label | Zero-based frames |
|---|---|---|
| 586875 | nogas, first run | excluded from the current analysis |
| 587178 | bc | 0, 5, 10, 15, 20, 25, 60, 110 |
| 587184 | spin | 0, 5, 10, 15, 20, 25, 55, 100 |
| 587207 | nogas, second run | 0, 1, 2, 3, 4, 5, 10, 25, 45 |
| 587214 | drop50 | 0, 25, 50, 57, 60, 71, 85, 99 (frame 60 added in September 2026) |
| 587217 | bcnogas1 | 0, 1, 2, 3, 4, 5, 6, 11, 25, 35 |
| 587225 | n2n2 | 0, 10, 19, 21, 24, 27, 30, 40, 70, 134 |
| 587241 | airn2 | 0, 9, 18, 20, 25, 30, 35, 40, 45, 50, 60, 80, 129 |
| 587250 | airn2ITO | 0, 11, 22, 30, 40, 75, 100, 129 |

Frame lists can be stored as suggestions in
`configs/secondary_samples/scan_registry.json`. A suggestion is not treated as
confirmed and cannot silently unlock fitting. Confirm frames either by filling
the registry's `key_frames` field or by passing `--frames` when building the
candidate configs.

Confirmation is particularly important for the late cooling in 587207 and
587217, the intermediate-temperature stages in 587241 and 587250, which nogas
run is primary versus replicate, and the correct wording for frames 0-50 of
587214.

## Stage 0: prepare and validate one scan (no fitting)

This stage reads the acquisition format described in the README's Scope
section. Copy `configs/secondary_samples/scan_registry.example.json` to
`scan_registry.json`, fill in the local source paths, then, for scan 587178:

```bash
python scripts/preparation/prepare_secondary_scan.py --scan 587178
```

This validates:

1. temperature-table and raw SRS frame counts;
2. processed frame 0 mapping to raw `frameNo` 1;
3. temperature alignment;
4. complete FR/IP/OOP frame columns;
5. identical q grids;
6. finite, positive d5i monitor values.

It writes raw and d5i-normalised fitting inputs, a frame manifest, a
provenance file and a key-frame suggestion CSV under
`local_results/prepared/<sample>/scan_<scan>/`. The frame manifest of scan
587178 is committed as `example_data/secondary_scan/scan587178_frame_manifest.csv`
so the later stages can be tried without the raw data.

The normalisation is a frame-wise multiplier:

`I_norm(frame,q) = I_raw(frame,q) * median(d5i) / d5i(frame)`

It affects fitted areas but not q position or FWHM within a frame.

## Stage 1: inspect the selected raw profiles without fitting

After confirming the frame list, plot the unsmoothed FR/IP/OOP profiles. For
587178:

```bash
python scripts/diagnostics/plot_key_frames.py \
  --processed-dir local_results/prepared/dimitar_bc/scan_587178 \
  --scan-id 587178 \
  --sample-label "Dimitar_bc" \
  --frames 0,5,10,15,20,25,60,110 \
  --outdir local_results/587178/key_frames
```

Use this figure to confirm candidate peaks and to decide whether a sample needs
an added, removed or split fitting component. Manual peak labels are starting
guesses, not final reported positions.

## Stage 2: build candidate configs

```bash
python scripts/preparation/build_secondary_sample_configs.py \
  --scan 587178 \
  --sample-label "Dimitar_bc" \
  --manifest local_results/prepared/dimitar_bc/scan_587178/scan587178_frame_manifest.csv \
  --frames 0,5,10,15,20,25,60,110 \
  --output-dir local_results/587178/candidate_configs
```

These configs deliberately inherit the limited candidate-window library used
for 587214 (`configs/secondary_samples/scan_587214/representative/`). They are
not assumed correct for another sample. Review and edit peak count, q bounds,
FWHM bounds and local background only where the raw key-frame plots justify a
change, and record what was changed: the committed 587178, 587184 and 587207
configurations carry a `candidate_revision` block for exactly this purpose.

## Stage 3: fit only the key frames

Run one cut at a time so a bad config is caught cheaply. FR example, using the
committed reviewed configuration:

```bash
python scripts/fitting/fit_giwaxs_series.py \
  local_results/prepared/dimitar_bc/scan_587178/scan587178_FR_111frames_norm.txt \
  --config configs/secondary_samples/scan_587178/scan587178_FR_candidate_keyframes.json \
  --outdir local_results/587178/key_frame_fits/FR \
  --stage selected \
  --frames 0,5,10,15,20,25,60,110 \
  --bootstrap 0
```

Repeat with IP and OOP after the FR output completes normally.

Review each `figures/fits_<window>.png` and require:

1. background follows the local baseline rather than swallowing broad signal;
2. total fit follows the visible data;
3. individual components correspond to visible peaks/shoulders;
4. residuals do not contain repeated peak-shaped structure;
5. centres and widths are not repeatedly at bounds;
6. neighbouring components do not swap identities;
7. proposed q shifts are supported by confidence intervals and adjacent frames.

If a window fails, revise that window and save a new versioned results folder.
Do not continue to all frames merely because the optimiser reports success.

## Stage 4: freeze the accepted profile family per window

After FR/IP/OOP key-frame fits are accepted, generate configs that lock each
window to its selected Gaussian, Lorentzian or pseudo-Voigt family. A selected
pseudo-Voigt eta is also locked.

```bash
python scripts/preparation/build_secondary_sample_configs.py \
  --scan 587178 \
  --sample-label "Dimitar_bc" \
  --manifest local_results/prepared/dimitar_bc/scan_587178/scan587178_frame_manifest.csv \
  --frames 0,5,10,15,20,25,60,110 \
  --freeze-from-results local_results/587178/key_frame_fits \
  --output-dir local_results/587178/locked_configs
```

This prevents the profile family from changing when the number of fitted frames
changes.

## Stage 5: decide whether all frames are necessary

Key-frame results are sufficient for samples used only as supporting condition
comparisons. Fit all frames only when the dissertation discusses the timing of
a transition, continuous q movement, or annealing kinetics.

If justified, run the locked config with `--bootstrap 0` first:

```bash
python scripts/fitting/fit_giwaxs_series.py \
  local_results/prepared/dimitar_bc/scan_587178/scan587178_FR_111frames_norm.txt \
  --config local_results/587178/locked_configs/scan587178_FR_locked_allframes.json \
  --outdir local_results/587178/all_frame_fits/FR \
  --stage all \
  --confirm-first-fit \
  --bootstrap 0
```

Run IP/OOP only where sector tracking materially supports the dissertation.
Do not perform 100-repeat bootstrap runs across every frame. Reserve bootstrap
for a small number of conclusions whose uncertainty matters.

## Stage 6: build conservative dissertation tables

For the key-frame results:

```bash
python scripts/reporting/build_secondary_sample_report.py \
  --scan 587178 \
  --sample-label "Dimitar_bc" \
  --manifest local_results/prepared/dimitar_bc/scan_587178/scan587178_frame_manifest.csv \
  --results-root local_results/587178/key_frame_fits \
  --outdir local_results/587178/key_frame_report
```

Use `accepted_peak_positions.csv` for the main q/d table. A final accepted
position must have a nonblank `q0_reported_Ainv`, a freely resolved
strong/detected status, delta BIC at least 6, area SNR at least 3, optimiser
success, and no centre/FWHM bound. `q0_fit_Ainv` is retained only for audit
and must not replace a rejected or blank reported position.

Fixed-shape fits may support peak presence/area when their detection gate
passes, but their borrowed q position must not be presented as an independent
measurement.

Components that fell just below the gate (delta BIC 2-6 or area SNR 2-3) are
listed in `tentative_or_excluded_peaks.csv`. Show them when a reviewer needs
to see that a visible feature was tested rather than ignored; nothing in that
table is an accepted position. Note that the fitter assigns the detection
status from the preliminary all-free fit and re-derives it only in the anchor
frame, so a tentative row can carry a final area SNR above 3.

Windows may carry `frame_min`/`frame_max` (zero-based, inclusive) in addition to
the temperature gate. Use this only for isothermal holds whose pattern changes
with time, as in the 587214 46 °C hold (frames 0-30 versus 31-52).

## Practical order

1. Use the completed 587178 key-frame report for peak assignment and the draft
   Results and Discussion chapter.
2. Do not expand 587178 to all 111 frames unless the dissertation needs the
   precise timing or kinetics of a transition.
3. 587214 visual QC is complete (v3, 9 September 2026); use that package.
4. Roll the settled key-frame workflow out to the remaining supporting samples
   one at a time.
5. Reserve profile freezing, all-frame fits and bootstrap uncertainty for a
   specific result that genuinely needs them.
