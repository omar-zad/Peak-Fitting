# Secondary-sample GIWAXS peak-fitting workflow

## Aim and scientific boundary

Use this workflow for the non-DROP40 samples to obtain reproducible peak
presence, fitted peak position (`q0`), d-spacing (`d = 2*pi/q0`), integrated
area, IP/OOP sector comparison, and changes with time/temperature.

Do not use this workflow to claim CCL, crystallite size, strain or disorder for
the secondary samples. FWHM remains a fitting diagnostic only.

No smoothing, vertical offset or MCMC is used. All initial runs use
`--bootstrap 0`.

## Status before new fitting

- Scan 586875 is intentionally excluded from the current dissertation
  peak-fitting workflow at the user's request.
- Scan 587214 already has prepared inputs and preliminary seven-key-frame fits.
  These remain a workflow test, not blanket scientific approval: their overlays
  and QC flags still require review.
- Scan 587178 has completed the eight confirmed key-frame screen in FR, IP and
  OOP using sample-specific candidate configs, `--bootstrap 0`, visual fit QC,
  and a conservative report. It is not an all-frame or CCL/FWHM analysis.
- Scan 587184 has completed the eight confirmed key-frame screen in FR, IP and
  OOP with sample-specific candidate-v2 configs, `--bootstrap 0`, full visual
  QC, an auditable manual exclusion, and a simple peaks-per-frame workbook in
  `results/dimitar_spin/scan_587184/key_frame_report_v2/`.
- Scan 587207 has completed the nine confirmed key-frame screen in FR, IP and
  OOP with temperature-regime-specific candidate-v5 configs, `--bootstrap 0`,
  full visual QC, auditable manual exclusions for unstable FR high-q component
  centres, and a simple peaks-per-frame workbook in
  `results/dimitar_nogas_587207/scan_587207/key_frame_report_v5/`. Frame 45 is
  labelled final/slight cooling, not a separately established cooled state.
- Scan 587217 has completed the ten confirmed key-frame screen in FR, IP and
  OOP with candidate-v2 configs and `--bootstrap 0`. One background-sensitive
  OOP shoulder was removed by visual QC, leaving 89 dissertation-facing q
  positions. Frame 35 is late cooling at 73.9 °C. Results and the simple
  workbook are in `results/dimitar_bcnogas1/scan_587217/key_frame_report_v2/`.
- Scan 587225 has completed the ten confirmed key-frame screen in FR, IP and
  OOP with revised candidate-v3 configs and `--bootstrap 0`. Reference-like
  underlays are excluded by construction and three dubious positions were
  removed by visual QC, leaving 135 dissertation-facing q positions. Results
  and the simple workbook are in
  `results/dimitar_n2n2/scan_587225/key_frame_report_v3/`.
- Scan 587241 has completed the thirteen confirmed key-frame screen in FR, IP
  and OOP with candidate-v1 configs and `--bootstrap 0`. The transition occurs
  between the selected 48.8 and 61.6 °C frames; 150 positions passed the
  conservative gate, with asymmetric or under-resolved lines retained only as
  cautious centre/presence observations. Results and the simple workbook are
  in `results/dimitar_airn2/scan_587241/key_frame_report_v1/`.
- Scan 587250 has completed the eight confirmed key-frame screen in FR, IP and
  OOP with candidate-v3 configs and `--bootstrap 0`. FR uses workbook sheet
  `pilatus2-587250_FULL_manual_wed`; the edited subset sheet remains excluded.
  Strong asymmetric high-q features are retained for provisional centre
  comparison only, not width interpretation. Results and the simple workbook
  are in `results/dimitar_airn2ito/scan_587250/key_frame_report_v3/`.

## Proposed key frames (review before fitting)

| Scan | Condition label | Proposed zero-based frames |
|---|---|---|
| 586875 | nogas, first run | excluded from the current analysis |
| 587178 | bc | 0, 5, 10, 15, 20, 25, 60, 110 (confirmed) |
| 587184 | spin | 0, 5, 10, 15, 20, 25, 55, 100 (confirmed) |
| 587207 | nogas, second run | 0, 1, 2, 3, 4, 5, 10, 25, 45 (confirmed) |
| 587214 | drop50 | 0, 25, 50, 57, 71, 85, 99 (already agreed) |
| 587217 | bcnogas1 | 0, 1, 2, 3, 4, 5, 6, 11, 25, 35 (confirmed) |
| 587225 | n2n2 | 0, 10, 19, 21, 24, 27, 30, 40, 70, 134 (confirmed) |
| 587241 | airn2 | 0, 9, 18, 20, 25, 30, 35, 40, 45, 50, 60, 80, 129 (confirmed) |
| 587250 | airn2ITO | 0, 11, 22, 30, 40, 75, 100, 129 (confirmed) |

The lists are stored as suggestions in
`configs/secondary_samples/scan_registry.json`. A suggestion is not treated as
confirmed and cannot silently unlock fitting. Confirm frames either by filling
the registry's `key_frames` field or by passing `--frames` when building the
candidate configs.

User/protocol confirmation is particularly important for the late cooling in
587207 and 587217, the intermediate-temperature stages in 587241 and 587250,
which nogas run is primary versus replicate, and the correct wording for frames
0-50 of 587214.

## Stage 0: prepare and validate one scan (no fitting)

Example for scan 587178:

```bash
python3 scripts/preparation/prepare_secondary_scan.py --scan 587178
```

This validates:

1. temperature-table and raw SRS frame counts;
2. processed frame 0 mapping to raw `frameNo` 1;
3. temperature alignment;
4. complete FR/IP/OOP frame columns;
5. identical q grids;
6. finite, positive d5i monitor values.

It writes raw and d5i-normalized fitting inputs, a manifest, source provenance,
and a key-frame suggestion CSV under:

`data/<sample>/scan_<scan>/processed/`

The normalization is a frame-wise multiplier:

`I_norm(frame,q) = I_raw(frame,q) * median(d5i) / d5i(frame)`

It affects fitted areas but not q position or FWHM within a frame.

## Stage 1: inspect the selected raw profiles without fitting

After confirming the frame list, plot the unsmoothed FR/IP/OOP profiles. For
587178:

```bash
python3 scripts/diagnostics/plot_key_frames.py \
  --processed-dir data/dimitar_bc/scan_587178/processed \
  --scan-id 587178 \
  --sample-label "Dimitar_bc" \
  --frames 0,5,10,15,20,25,60,110 \
  --outdir results/dimitar_bc/scan_587178/key_frames_v2
```

Use this figure to confirm candidate peaks and to decide whether a sample needs
an added, removed or split fitting component. Manual peak labels are starting
guesses, not final reported positions.

## Stage 2: build candidate configs

```bash
python3 scripts/preparation/build_secondary_sample_configs.py \
  --scan 587178 \
  --sample-label "Dimitar_bc" \
  --manifest data/dimitar_bc/scan_587178/processed/scan587178_frame_manifest.csv \
  --frames 0,5,10,15,20,25,60,110 \
  --output-dir configs/secondary_samples/scan_587178/candidate_v2
```

These configs deliberately inherit the limited candidate-window library used
for 587214. They are not assumed correct for another sample. Review and edit
peak count, q bounds, FWHM bounds and local background only where the raw
key-frame plots justify a change.

## Stage 3: fit only the key frames

Run one cut at a time so a bad config is caught cheaply. FR example:

```bash
python3 scripts/fitting/fit_giwaxs_series.py \
  data/dimitar_bc/scan_587178/processed/scan587178_FR_111frames_norm.txt \
  --config configs/secondary_samples/scan_587178/candidate_v2/scan587178_FR_candidate_keyframes.json \
  --outdir results/dimitar_bc/scan_587178/key_frame_fits_v2/FR \
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
python3 scripts/preparation/build_secondary_sample_configs.py \
  --scan 587178 \
  --sample-label "Dimitar_bc" \
  --manifest data/dimitar_bc/scan_587178/processed/scan587178_frame_manifest.csv \
  --frames 0,5,10,15,20,25,60,110 \
  --freeze-from-results results/dimitar_bc/scan_587178/key_frame_fits_v2 \
  --output-dir configs/secondary_samples/scan_587178/locked
```

This prevents the profile family from changing when the number of fitted frames
changes.

## Stage 5: decide whether all frames are necessary

Key-frame results are sufficient for samples used only as supporting condition
comparisons. Fit all frames only when the dissertation discusses the timing of
a transition, continuous q movement, or annealing kinetics.

If justified, run the locked config with `--bootstrap 0` first:

```bash
python3 scripts/fitting/fit_giwaxs_series.py \
  data/dimitar_bc/scan_587178/processed/scan587178_FR_111frames_norm.txt \
  --config configs/secondary_samples/scan_587178/locked/scan587178_FR_locked_allframes.json \
  --outdir results/dimitar_bc/scan_587178/all_frame_fits_v1/FR \
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
python3 scripts/reporting/build_secondary_sample_report.py \
  --scan 587178 \
  --sample-label "Dimitar_bc" \
  --manifest data/dimitar_bc/scan_587178/processed/scan587178_frame_manifest.csv \
  --results-root results/dimitar_bc/scan_587178/key_frame_fits_v2 \
  --outdir results/dimitar_bc/scan_587178/key_frame_report_v2
```

Use `accepted_peak_positions.csv` for the main q/d table. A final accepted
position must have a nonblank `q0_reported_Ainv`, a freely resolved
strong/detected status, ΔBIC at least 6, area SNR at least 3, optimiser success,
and no centre/FWHM bound. `q0_fit_Ainv` is retained only for audit and must not
replace a rejected or blank reported position.

Fixed-shape fits may support peak presence/area when their detection gate
passes, but their borrowed q position must not be presented as an independent
measurement.

## Practical order

1. Use the completed 587178 key-frame report for peak assignment and the draft
   Results and Discussion chapter.
2. Do not expand 587178 to all 111 frames unless the dissertation needs the
   precise timing or kinetics of a transition.
3. Finish visual QC of the existing 587214 key-frame fits.
4. Roll the settled key-frame workflow out to the remaining supporting samples
   one at a time.
5. Reserve profile freezing, all-frame fits and bootstrap uncertainty for a
   specific result that genuinely needs them.
