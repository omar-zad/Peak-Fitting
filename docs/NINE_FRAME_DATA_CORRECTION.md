# Drop40 nine-frame data correction

The final Drop40 temperature series is:

| Frame | Full scan | Nominal temperature (°C) | Measured `temp3` (°C) |
|---:|---:|---:|---:|
| 193 | 587193 | 40 | 39.9 |
| 194 | 587194 | 50 | 49.9 |
| 195 | 587195 | 55 | 54.8 |
| 196 | 587196 | 60 | 59.9 |
| 197 | 587197 | 65 | 64.5 |
| 198 | 587198 | 70 | 69.6 |
| 199 | 587199 | 75 | 74.5 |
| 200 | 587200 | 80 | 79.2 |
| 201 | 587201 | 150 | 145.4 |

The old `data/drop40/provenance/drop40_FR.txt` contains frames 193–200 only. Its previous fitting
configuration incorrectly associated frame 199 with 80 °C and frame 200 with
150 °C. It is retained solely for provenance.

Use the corrected inputs in `data/drop40/processed/`:

- `drop40_FR_9frames.txt`
- `drop40_IP_9frames.txt`
- `drop40_OOP_9frames.txt`
- `drop40_FR_9frames_norm.txt`
- `drop40_IP_9frames_norm.txt`
- `drop40_OOP_9frames_norm.txt`
- `drop40_9frame_manifest.csv`

Regenerate them from the original sector CSVs with:

```bash
python3 scripts/preparation/build_drop40_nine_frame_inputs.py
```

The source exposure time and transmission are constant across the nine scans.
The norm-corrected files multiply each scan by
`median(norm across all scans) / scan norm`. This preserves approximately the
original intensity scale while correcting the roughly 6% recorded monitor
variation. Raw and normalized results should be compared before interpreting
areas quantitatively.

All peak-fitting outputs, spreadsheets, and slides created before this
correction are superseded and must not be used as final numerical results.
